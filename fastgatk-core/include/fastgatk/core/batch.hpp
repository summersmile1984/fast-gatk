#pragma once

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace fastgatk::core {

// HostBatch is deliberately metadata-only: individual kernels own their strongly
// typed std::vector/arena payloads, while this common object gives the runtime one
// accounting and schema boundary for all of them.
struct HostBatch {
    std::string schema;
    std::size_t records = 0;
    std::size_t bytes = 0;

    explicit HostBatch(std::string schema_name = {}) : schema(std::move(schema_name)) {}

    template<class Container>
    void observe(const Container& values, std::size_t logical_records = 0) {
        using Value = typename Container::value_type;
        bytes += values.size() * sizeof(Value);
        if (logical_records != 0) records = logical_records;
    }
};

template<class ExecSpace>
class DeviceBatch {
public:
    using execution_space = ExecSpace;
    using memory_space = typename ExecSpace::memory_space;

    explicit DeviceBatch(std::size_t records = 0) : records_(records) {}

    template<class View>
    void bind(std::string_view label, const View& view) {
        bind(label, view, view.span());
    }

    // Bind a reusable capacity view while accounting only the logical prefix
    // consumed by this batch.  Persistent plans can therefore retain device
    // allocations between batches without overstating per-batch telemetry.
    template<class View>
    void bind(std::string_view label, const View& view, std::size_t logical_span) {
        using ViewType = std::remove_cv_t<std::remove_reference_t<View>>;
        static_assert(std::is_same_v<typename ViewType::memory_space, memory_space>,
                      "DeviceBatch View must belong to the plan execution memory space");
        bindings_.emplace_back(label);
        const auto span = logical_span < view.span() ? logical_span : view.span();
        bytes_ += span * sizeof(typename ViewType::value_type);
    }

    std::size_t records() const { return records_; }
    std::size_t bytes() const { return bytes_; }
    std::size_t bindings() const { return bindings_.size(); }
    const std::vector<std::string>& labels() const { return bindings_; }

private:
    std::size_t records_ = 0;
    std::size_t bytes_ = 0;
    std::vector<std::string> bindings_;
};

}  // namespace fastgatk::core
