#pragma once

#define LOG_TRACE(...) spdlog::get("stderr")->trace(__VA_ARGS__)        // NOLINT
#define LOG_DEBUG(...) spdlog::get("stderr")->debug(__VA_ARGS__)        // NOLINT
#define LOG_INFO(...) spdlog::get("stderr")->info(__VA_ARGS__)          // NOLINT
#define LOG_WARN(...) spdlog::get("stderr")->warn(__VA_ARGS__)          // NOLINT
#define LOG_ERROR(...) spdlog::get("stderr")->error(__VA_ARGS__)        // NOLINT
#define LOG_CRITICAL(...) spdlog::get("stderr")->critical(__VA_ARGS__)  // NOLINT
