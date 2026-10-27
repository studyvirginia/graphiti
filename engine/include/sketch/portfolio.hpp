#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <memory>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#ifdef __APPLE__
#include <pthread/qos.h>
#endif

#include "sketch/multitree_gp.hpp"

namespace Sketch {

struct PortfolioMemberResult {
    MultitreeResult result{};
    double seconds{0.0};
    bool ok{false};
    std::string error;
};

#ifdef __APPLE__
#include <sys/sysctl.h>
inline auto getPerformanceCoreCount() -> int {
    int pcores = 0;
    size_t len = sizeof(pcores);
    if (sysctlbyname("hw.perflevel0.logicalcpu", &pcores, &len, nullptr, 0) == 0 && pcores > 0) {
        return pcores;
    }
    if (sf2On(sf2Mask(), 14)) { return std::max(1, static_cast<int>(std::thread::hardware_concurrency())); }
    return 2;
}
#else
inline auto getPerformanceCoreCount() -> int {
    return std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
}
#endif

template <typename DTable>
std::vector<PortfolioMemberResult> runPortfolio(
    std::vector<Operon::Scalar> const& dx,
    std::vector<Operon::Scalar> const& dy,
    std::vector<MultitreeConfig> const& members,
    NoiseModel const& noise,
    int nThreads,
    bool cppHalving = false,
    std::size_t scoutGens = 10,
    std::size_t midGens = 25,
    int nSurvivors = 0)
{
    std::vector<PortfolioMemberResult> out(members.size());
    if (members.empty()) { return out; }

    auto const pCores = getPerformanceCoreCount();
    auto const effectiveThreads = (nThreads <= 0) ? pCores : std::min(nThreads, pCores);
    auto const nWorkers = static_cast<std::size_t>(
        std::max(1, std::min<int>(effectiveThreads, static_cast<int>(members.size()))));

    if (!cppHalving || members.size() <= 1) {
        std::atomic<std::size_t> next{0};
        auto worker = [&]() {
#ifdef __APPLE__
            pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
            DTable dtable;
            for (;;) {
                auto const i = next.fetch_add(1, std::memory_order_relaxed);
                if (i >= members.size()) { break; }
                auto const t0 = std::chrono::steady_clock::now();
                try {
                    MultitreeGP<DTable> gp(dtable, dx, dy, members[i], noise);
                    out[i].result = gp.Run();
                    out[i].ok = true;
                } catch (std::exception const& e) {
                    out[i].ok = false;
                    out[i].error = e.what();
                } catch (...) {
                    out[i].ok = false;
                    out[i].error = "unknown exception";
                }
                out[i].seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t0).count();
            }
        };

        if (nWorkers == 1) {
            worker();
            return out;
        }

        std::vector<std::thread> pool;
        pool.reserve(nWorkers);
        for (std::size_t w = 0; w < nWorkers; ++w) { pool.emplace_back(worker); }
        for (auto& th : pool) { th.join(); }
        return out;
    }

    std::vector<std::unique_ptr<MultitreeGP<DTable>>> engines;
    engines.reserve(members.size());
    DTable dummyDTable;
    for (std::size_t i = 0; i < members.size(); ++i) {
        engines.emplace_back(std::make_unique<MultitreeGP<DTable>>(
            dummyDTable, dx, dy, members[i], noise));
    }

    auto runParallelTasks = [&](std::size_t numTasks, auto const& taskFn) {
        if (nWorkers <= 1 || numTasks <= 1) {
            DTable dt;
            for (std::size_t i = 0; i < numTasks; ++i) {
                taskFn(i, dt);
            }
            return;
        }
        std::atomic<std::size_t> nextTask{0};
        auto worker = [&]() {
#ifdef __APPLE__
            pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
            DTable dt;
            for (;;) {
                auto const taskIdx = nextTask.fetch_add(1, std::memory_order_relaxed);
                if (taskIdx >= numTasks) { break; }
                taskFn(taskIdx, dt);
            }
        };
        std::vector<std::thread> pool;
        pool.reserve(nWorkers);
        for (std::size_t w = 0; w < nWorkers; ++w) { pool.emplace_back(worker); }
        for (auto& th : pool) { th.join(); }
    };

    std::vector<char> memberOk(members.size(), 0);
    std::vector<std::string> memberError(members.size(), "");
    std::vector<double> memberSeconds(members.size(), 0.0);

    bool const useTriWave = (midGens > scoutGens) && (members.size() > 2);
    std::size_t K1 = 0;
    std::size_t K2 = 0;
    if (useTriWave) {
        if (nSurvivors > 0) {
            K2 = std::min<std::size_t>(members.size(), static_cast<std::size_t>(nSurvivors));
            K1 = std::max<std::size_t>(K2, (members.size() + 1) / 2);
        } else {
            K1 = std::max<std::size_t>(1, (members.size() + 1) / 2);
            K2 = std::max<std::size_t>(1, (K1 + 1) / 2);
        }
    } else {
        K1 = (nSurvivors > 0)
            ? std::min<std::size_t>(members.size(), static_cast<std::size_t>(nSurvivors))
            : std::max<std::size_t>(1, members.size() / 2);
    }

    std::vector<std::size_t> active(members.size());
    std::iota(active.begin(), active.end(), 0);
    std::vector<std::size_t> pruned;

    runParallelTasks(active.size(), [&](std::size_t taskIdx, DTable& dt) {
        auto const i = active[taskIdx];
        auto const t0 = std::chrono::steady_clock::now();
        try {
            engines[i]->SetDispatchTable(&dt);
            engines[i]->InitSearchState();
            engines[i]->RunTo(std::min(scoutGens, members[i].generations));
            memberOk[i] = true;
        } catch (std::exception const& e) {
            memberOk[i] = false;
            memberError[i] = e.what();
        } catch (...) {
            memberOk[i] = false;
            memberError[i] = "unknown exception in scout phase";
        }
        memberSeconds[i] += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    });

    std::stable_sort(active.begin(), active.end(), [&](std::size_t a, std::size_t b) {
        if (memberOk[a] != memberOk[b]) { return memberOk[a] > memberOk[b]; }
        return engines[a]->BestShapeError() < engines[b]->BestShapeError();
    });

    std::size_t const k1 = std::min(K1, active.size());
    for (std::size_t idx = k1; idx < active.size(); ++idx) {
        pruned.push_back(active[idx]);
    }
    active.resize(k1);

    if (useTriWave && active.size() > 1) {
        runParallelTasks(active.size(), [&](std::size_t taskIdx, DTable& dt) {
            auto const i = active[taskIdx];
            if (!memberOk[i]) { return; }
            auto const t0 = std::chrono::steady_clock::now();
            try {
                engines[i]->SetDispatchTable(&dt);
                engines[i]->RunTo(std::min(midGens, members[i].generations));
                memberOk[i] = true;
            } catch (std::exception const& e) {
                memberOk[i] = false;
                memberError[i] = e.what();
            } catch (...) {
                memberOk[i] = false;
                memberError[i] = "unknown exception in mid wave";
            }
            memberSeconds[i] += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
        });

        std::stable_sort(active.begin(), active.end(), [&](std::size_t a, std::size_t b) {
            if (memberOk[a] != memberOk[b]) { return memberOk[a] > memberOk[b]; }
            return engines[a]->BestShapeError() < engines[b]->BestShapeError();
        });

        std::size_t const k2 = std::min(K2, active.size());
        for (std::size_t idx = k2; idx < active.size(); ++idx) {
            pruned.push_back(active[idx]);
        }
        active.resize(k2);
    }

    std::size_t const survivorThreads = sf2On(sf2Mask(), 14)
        ? (active.empty() ? std::size_t{1} : std::max<std::size_t>(1, static_cast<std::size_t>(pCores) / active.size()))
        : ((pCores >= 4 && active.size() > 0) ? static_cast<std::size_t>(pCores / active.size()) : 1);
    for (auto i : active) {
        if (survivorThreads > 1) {
            engines[i]->SetThreads(survivorThreads);
        }
    }

    runParallelTasks(active.size(), [&](std::size_t taskIdx, DTable& dt) {
        auto const i = active[taskIdx];
        auto const t0 = std::chrono::steady_clock::now();
        if (!memberOk[i]) {
            out[i].ok = false;
            out[i].error = memberError[i];
        } else {
            try {
                engines[i]->SetDispatchTable(&dt);
                engines[i]->RunTo(members[i].generations);
                out[i].result = engines[i]->Finalize();
                out[i].ok = true;
            } catch (std::exception const& e) {
                out[i].ok = false;
                out[i].error = e.what();
            } catch (...) {
                out[i].ok = false;
                out[i].error = "unknown exception in final wave";
            }
        }
        memberSeconds[i] += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    });

    runParallelTasks(pruned.size(), [&](std::size_t taskIdx, DTable& dt) {
        auto const i = pruned[taskIdx];
        auto const t0 = std::chrono::steady_clock::now();
        if (!memberOk[i]) {
            out[i].ok = false;
            out[i].error = memberError[i];
        } else {
            try {
                engines[i]->SetDispatchTable(&dt);
                out[i].result = engines[i]->Finalize();
                out[i].ok = true;
            } catch (std::exception const& e) {
                out[i].ok = false;
                out[i].error = e.what();
            } catch (...) {
                out[i].ok = false;
                out[i].error = "unknown exception in finalize";
            }
        }
        memberSeconds[i] += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    });

    for (std::size_t i = 0; i < members.size(); ++i) {
        out[i].seconds = memberSeconds[i];
    }

    return out;
}

}
