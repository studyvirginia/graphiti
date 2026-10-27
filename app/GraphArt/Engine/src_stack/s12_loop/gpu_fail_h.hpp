#pragma once
#include <stdexcept>
[[noreturn]] inline void gpu_fail() { throw std::runtime_error("Metal command buffer failed"); }
