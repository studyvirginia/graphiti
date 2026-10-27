#ifndef SKETCH_MULTITREE_GP_HPP
#define SKETCH_MULTITREE_GP_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <functional>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <thread>
#include <vector>
#include <unordered_set>

#ifdef SKETCH_DUP_STATS
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#endif

#include <operon/operon/core/dataset.hpp>
#include <operon/operon/core/node.hpp>
#include <operon/operon/core/problem.hpp>
#include <operon/operon/core/pset.hpp>
#ifdef __APPLE__
#include <pthread.h>
#include <pthread/qos.h>
#include <Accelerate/Accelerate.h>
#endif
#include <operon/operon/interpreter/interpreter.hpp>
#include <operon/operon/operators/creator.hpp>
#include <operon/operon/operators/crossover.hpp>
#include <operon/operon/operators/initializer.hpp>
#include <operon/operon/operators/mutation.hpp>
#include <operon/operon/optimizer/optimizer.hpp>
#include <operon/operon/formatter/formatter.hpp>

#include "pair_individual.hpp"
#include "grid_sdf.hpp"
#include "bytecode_evaluator.hpp"
#include "register_interp.hpp"
#include "shape_const_optimizer.hpp"
#include "point_segment_kernel.hpp"
#include "one_residual.hpp"
#include "persistent_worker_pool.hpp"
#include "batched_pointseg_gpu.hpp"
#include "concurrent_lm_gpu.hpp"
#include "cis_scoring.hpp"
#include "sigma_v2_native.hpp"
#include "symmetric_fitness.hpp"
#include "robust_continuous_scoring.hpp"
#include "frequency_discovery.hpp"
#include "matrix_pencil.hpp"
#include "sf2_guard.hpp"
#include "op2_guard.hpp"
#include "topology_detect.hpp"
#include "bracket_newton_runtime.hpp"
#include "pair_amplitude_adapter.hpp"
#include "stroke_downsample.hpp"
#include "parametrization_gauge.hpp"
#include "chow_abstention.hpp"
#include "nml_parsimony.hpp"
#include "t13_guard.hpp"

#include "multitree_gp_types.hpp"

template <typename DTable>
class MultitreeGP {
public:
#include "multitree_gp_construct.hpp"
#include "multitree_gp_eval.hpp"
#include "multitree_gp_finalize.hpp"
#include "multitree_gp_scoring.hpp"
#include "multitree_gp_cis.hpp"
#include "multitree_gp_driver.hpp"
};

}
#endif
