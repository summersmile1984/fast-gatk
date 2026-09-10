#include "fastgatk/runtime/resource.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <cstdlib>

int main() {
    try {
        const auto snapshot = fastgatk::runtime::ResourceSnapshot::probe();
        if (snapshot.hardware_threads == 0 || snapshot.effective_threads() == 0 ||
            snapshot.effective_threads(1) != 1)
            throw std::runtime_error("invalid effective resource limits");
        const auto json = snapshot.to_json();
        if (json.find("hardware_threads") == std::string::npos ||
            json.find("memory_limit_bytes") == std::string::npos ||
            json.find("scratch_free_bytes") == std::string::npos ||
            json.find("scratch_free_inodes") == std::string::npos ||
            json.find("file_descriptor_limit") == std::string::npos ||
            json.find("host_hard_bytes") == std::string::npos ||
            json.find("device_memory_bytes") == std::string::npos ||
            json.find("device_free_bytes") == std::string::npos ||
            json.find("reported_device_memory_bytes") == std::string::npos ||
            json.find("device_backend_compiled") == std::string::npos ||
            json.find("device_assignment_available") == std::string::npos ||
            json.find("device_telemetry_rejected") == std::string::npos ||
            json.find("device_runtime_probe") == std::string::npos ||
            json.find("scratch_hard_bytes") == std::string::npos ||
            json.find("cpu_threads") == std::string::npos ||
            json.find("local_ssd") == std::string::npos ||
            json.find("remote_input") == std::string::npos)
            throw std::runtime_error("resource JSON schema missing required fields");
        if (!snapshot.memory_budget_allows(0) ||
            (snapshot.safe_memory_budget_bytes() > snapshot.effective_memory_limit_bytes() &&
             snapshot.effective_memory_limit_bytes() != 0))
            throw std::runtime_error("invalid safe memory budget");

        const fastgatk::runtime::ResourceBudget budget{
            1000, 1000, 1000, 100, 100};
        const fastgatk::runtime::WorkEstimate work{
            64, 8, 100, 200, 100};
        const fastgatk::runtime::AdaptiveController controller;
        const auto initial = controller.initial(budget, work);
        if (initial.max_reads != 3 || initial.max_host_bytes != 600 ||
            initial.max_device_bytes != 600)
            throw std::runtime_error("adaptive initial limits are not conservative");
        fastgatk::runtime::RuntimeTelemetry high_memory;
        high_memory.host_bytes = 500;
        auto reduced = controller.next(high_memory, initial,
                                       fastgatk::runtime::Pressure::HostMemory);
        if (reduced.max_reads != 1)
            throw std::runtime_error("adaptive controller did not reduce under pressure");
        fastgatk::runtime::RuntimeTelemetry low_memory;
        low_memory.compute_queue_empty = true;
        const auto grown = controller.next(low_memory, reduced,
                                           fastgatk::runtime::Pressure::Normal);
        if (grown.max_reads <= reduced.max_reads)
            throw std::runtime_error("adaptive controller did not grow an idle batch");
        fastgatk::runtime::RuntimeTelemetry scratch_pressure;
        scratch_pressure.scratch_bytes = 801;
        if (!controller.should_spill(scratch_pressure, budget) ||
            !controller.should_pause_reader(scratch_pressure, budget))
            throw std::runtime_error("scratch pressure was not propagated");
        if (budget.host_hard_bytes != 0 || budget.device_hard_bytes != 0 ||
            budget.scratch_hard_bytes != 0 || budget.cpu_threads != 0 ||
            budget.local_ssd || budget.remote_input)
            throw std::runtime_error("legacy budget aggregate compatibility failed");

#if defined(__linux__)
        // A scheduler can expose CUDA selectors to a CPU-only tool process.
        // Preserve that provenance but never admit its memory budget unless a
        // device execution space was compiled into this binary.
        if (setenv("CUDA_VISIBLE_DEVICES", "0,1", 1) != 0 ||
            setenv("FASTGATK_DEVICE_MEMORY_BYTES", "4096", 1) != 0 ||
            setenv("FASTGATK_DEVICE_FREE_BYTES", "2048", 1) != 0)
            throw std::runtime_error("unable to configure synthetic device allocation");
        const auto device_snapshot = fastgatk::runtime::ResourceSnapshot::probe();
        if (device_snapshot.visible_gpus != 2 ||
            device_snapshot.device_visibility_source != "CUDA_VISIBLE_DEVICES" ||
            device_snapshot.reported_device_memory_bytes != 4096 ||
            device_snapshot.reported_device_free_bytes != 2048)
            throw std::runtime_error("device allocation provenance probe failed");
        if (!device_snapshot.device_backend_compiled) {
            if (device_snapshot.device_assignment_available ||
                !device_snapshot.device_telemetry_rejected ||
                device_snapshot.device_memory_bytes != 0 ||
                device_snapshot.device_free_bytes != 0)
                throw std::runtime_error("host build admitted device budget");
        } else if (!device_snapshot.device_assignment_available ||
                   device_snapshot.device_telemetry_rejected ||
                   device_snapshot.device_memory_bytes != 4096 ||
                   device_snapshot.device_free_bytes != 2048) {
            throw std::runtime_error("device build rejected assigned device budget");
        }
#endif

        fastgatk::runtime::BoundedByteQueue<int> queue(10);
        if (!queue.push(7, 6) || queue.queued_bytes() != 6)
            throw std::runtime_error("byte queue push accounting failed");
        std::uint64_t popped_bytes = 0;
        const auto item = queue.pop(&popped_bytes);
        if (!item || *item != 7 || popped_bytes != 6 || queue.queued_bytes() != 0)
            throw std::runtime_error("byte queue pop accounting failed");
        queue.close();
        if (queue.push(8, 1) || queue.pop())
            throw std::runtime_error("closed byte queue accepted work");
        std::cout << "{\"status\":\"pass\",\"hardware_threads\":"
                  << snapshot.hardware_threads << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
