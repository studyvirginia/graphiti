#pragma once


#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

#ifndef _WIN32
#include <pthread.h>
#endif

namespace Sketch {

inline auto persistentPoolEnabledFromEnv() -> bool
{
    static bool const v = [] {
        char const* e = std::getenv("SKETCH_POOL");
        return (e != nullptr) && (std::atoi(e) != 0);
    }();
    return v;
}

#ifndef _WIN32
class PersistentWorkerPool {
public:
    static auto instance() -> PersistentWorkerPool&
    {
        static PersistentWorkerPool pool;
        return pool;
    }

    auto tryAcquireDispatcher() -> std::unique_lock<std::mutex>
    {
        return std::unique_lock<std::mutex>(dispatchMutex_, std::try_to_lock);
    }

    void dispatch(std::size_t workers, std::function<void()> const& body)
    {
        ensureCapacity(workers);
        workers = std::min(workers, capacity_);
        if (workers == 0) {
            body();
            return;
        }
        std::unique_lock<std::mutex> lk(m_);
        job_ = &body;
        activeWorkers_ = workers;
        remaining_ = workers;
        ++epoch_;
        cv_.notify_all();
        doneCv_.wait(lk, [this] { return remaining_ == 0; });
        job_ = nullptr;
    }

    ~PersistentWorkerPool()
    {
        {
            std::lock_guard<std::mutex> lk(m_);
            shutdown_ = true;
        }
        cv_.notify_all();
        for (pthread_t t : threads_) { pthread_join(t, nullptr); }
    }

    PersistentWorkerPool(PersistentWorkerPool const&) = delete;
    PersistentWorkerPool& operator=(PersistentWorkerPool const&) = delete;

private:
    PersistentWorkerPool() = default;

    struct WorkerCtx {
        PersistentWorkerPool* pool;
        std::size_t idx;
    };

    static void* workerTrampoline(void* p)
    {
        std::unique_ptr<WorkerCtx> ctx(static_cast<WorkerCtx*>(p));
        ctx->pool->workerLoop(ctx->idx);
        return nullptr;
    }

    void workerLoop(std::size_t idx)
    {
#ifdef __APPLE__
        pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
        std::uint64_t seenEpoch = 0;
        for (;;) {
            std::function<void()> const* jobCopy = nullptr;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [this, idx, seenEpoch] {
                    return shutdown_ || (epoch_ != seenEpoch && idx < activeWorkers_);
                });
                if (shutdown_) { return; }
                seenEpoch = epoch_;
                jobCopy = job_;
            }
            (*jobCopy)();
            {
                std::lock_guard<std::mutex> lk(m_);
                --remaining_;
                if (remaining_ == 0) { doneCv_.notify_all(); }
            }
        }
    }

    void ensureCapacity(std::size_t n)
    {
        constexpr std::size_t kWorkerStackBytes = std::size_t{16} << 20;
        while (capacity_ < n) {
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setstacksize(&attr, kWorkerStackBytes);
            auto* ctx = new WorkerCtx{this, capacity_};
            pthread_t tid{};
            if (pthread_create(&tid, &attr, &PersistentWorkerPool::workerTrampoline, ctx) == 0) {
                threads_.push_back(tid);
                ++capacity_;
            } else {
                delete ctx;
                break;
            }
            pthread_attr_destroy(&attr);
        }
    }

    std::mutex dispatchMutex_;
    std::mutex m_;
    std::condition_variable cv_;
    std::condition_variable doneCv_;
    std::vector<pthread_t> threads_;
    std::size_t capacity_ = 0;
    std::size_t activeWorkers_ = 0;
    std::uint64_t epoch_ = 0;
    std::size_t remaining_ = 0;
    std::function<void()> const* job_ = nullptr;
    bool shutdown_ = false;
};
#endif

}
