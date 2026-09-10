#pragma once

#include "fastgatk/core/batch.hpp"

#include <Kokkos_Core.hpp>

#include <chrono>
#include <cstddef>
#include <string>
#include <string_view>

namespace fastgatk::core {

struct PlanTelemetry {
    std::string label;
    std::size_t prepare_bytes = 0;
    std::size_t device_bytes = 0;
    std::size_t device_bindings = 0;
    double prepare_seconds = 0.0;
    double execute_seconds = 0.0;
    std::size_t execute_calls = 0;
};

// KernelPlan is the common lifetime/accounting object. Module plans keep their
// typed Kokkos::View members, and use this class to make preparation and repeated
// execution explicit. This is what lets the runtime amortize packing/deep_copy and
// report kernel-only versus end-to-end performance consistently.
template<class ExecSpace>
class KernelPlan {
public:
    using execution_space = ExecSpace;
    using memory_space = typename ExecSpace::memory_space;

    explicit KernelPlan(std::string label) : telemetry_{std::move(label)} {}

    void begin_prepare(const HostBatch& host) {
        telemetry_.prepare_bytes = host.bytes;
        prepare_start_ = std::chrono::steady_clock::now();
    }

    void end_prepare(const DeviceBatch<ExecSpace>& device) {
        telemetry_.device_bytes = device.bytes();
        telemetry_.device_bindings = device.bindings();
        telemetry_.prepare_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - prepare_start_).count();
    }

    void begin_execute() { execute_start_ = std::chrono::steady_clock::now(); }

    void end_execute() {
        telemetry_.execute_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - execute_start_).count();
        ++telemetry_.execute_calls;
    }

    const PlanTelemetry& telemetry() const { return telemetry_; }
    const std::string& label() const { return telemetry_.label; }

private:
    PlanTelemetry telemetry_;
    std::chrono::steady_clock::time_point prepare_start_{};
    std::chrono::steady_clock::time_point execute_start_{};
};

}  // namespace fastgatk::core
