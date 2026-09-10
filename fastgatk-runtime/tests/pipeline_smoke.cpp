#include "fastgatk/runtime/pipeline.hpp"

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
    using Pipeline = fastgatk::runtime::ThreeStagePipeline<int, int, std::string>;
    int next = 0;
    std::vector<std::string> output;
    Pipeline pipeline(
        Pipeline::Limits{8, 8, 16},
        [&]() -> std::optional<int> {
            if (next == 64) return std::nullopt;
            return next++;
        },
        [](int value) -> std::optional<int> { return value * 3; },
        [](int value) -> std::optional<std::string> {
            return std::to_string(value);
        },
        [&](std::string value) { output.push_back(std::move(value)); },
        [](const int&) { return std::uint64_t{1}; },
        [](const int&) { return std::uint64_t{1}; },
        [](const std::string& value) {
            return static_cast<std::uint64_t>(value.size());
        });
    const auto metrics = pipeline.run();
    if (metrics.decoded_items != 64 || metrics.computed_items != 64 ||
        metrics.encoded_items != 64 || output.size() != 64 ||
        metrics.decoded_bytes != 64 || metrics.computed_bytes != 64 ||
        metrics.encoded_bytes == 0 || metrics.peak_decoded_bytes == 0 ||
        metrics.peak_computed_bytes == 0 || metrics.peak_encoded_bytes == 0)
        throw std::runtime_error("three-stage pipeline did not preserve all items");
    for (std::size_t index = 0; index < output.size(); ++index) {
        if (output[index] != std::to_string(static_cast<int>(index) * 3))
            throw std::runtime_error("three-stage pipeline changed output order");
    }

    bool rejected = false;
    try {
        Pipeline invalid(
            Pipeline::Limits{0, 1, 1},
            []() -> std::optional<int> { return std::nullopt; },
            [](int value) -> std::optional<int> { return value; },
            [](int value) -> std::optional<std::string> { return std::to_string(value); },
            [](std::string) {},
            [](const int&) { return std::uint64_t{1}; },
            [](const int&) { return std::uint64_t{1}; },
            [](const std::string& value) { return static_cast<std::uint64_t>(value.size()); });
        (void)invalid;
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    if (!rejected) throw std::runtime_error("zero-capacity pipeline was accepted");

    int failure_next = 0;
    bool propagated = false;
    try {
        Pipeline failing(
            Pipeline::Limits{8, 8, 8},
            [&]() -> std::optional<int> {
                if (failure_next == 32) return std::nullopt;
                return failure_next++;
            },
            [](int value) -> std::optional<int> {
                if (value == 7) throw std::runtime_error("synthetic compute failure");
                return value;
            },
            [](int value) -> std::optional<std::string> { return std::to_string(value); },
            [](std::string) {},
            [](const int&) { return std::uint64_t{1}; },
            [](const int&) { return std::uint64_t{1}; },
            [](const std::string& value) { return static_cast<std::uint64_t>(value.size()); });
        (void)failing.run();
    } catch (const std::runtime_error& error) {
        propagated = std::string(error.what()).find("synthetic compute failure") !=
                     std::string::npos;
    }
    if (!propagated) throw std::runtime_error("pipeline stage exception was not propagated");

    std::cout << "{\"status\":\"pass\",\"decoded\":" << metrics.decoded_items
              << ",\"peak_decoded_bytes\":" << metrics.peak_decoded_bytes << "}\n";
    return 0;
}
