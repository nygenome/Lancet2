#ifndef SRC_LANCET_BASE_LOGGING_H_
#define SRC_LANCET_BASE_LOGGING_H_

#include <type_traits>
#include <utility>

#include "lancet/base/types.h"
#include "spdlog/async.h"
#include "spdlog/async_logger.h"
#include "spdlog/common.h"
#include "spdlog/logger.h"
#include "spdlog/sinks/sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

namespace lancet {

static constexpr auto LOGGER_NAME = "LANCET_LOG";

template <typename Logger = spdlog::async_logger, typename Sink = spdlog::sinks::stderr_color_sink_mt, class... Args>
void RegisterLancetLogger(Args&&... args) {
  static_assert(std::is_base_of_v<spdlog::logger, Logger>, "Logger must have base of spdlog::logger");
  static_assert(std::is_base_of_v<spdlog::sinks::sink, Sink>, "Sink must implement spdlog::sinks::sink interface");

  // Already registered logger previously, so we can return early
  // NOLINTNEXTLINE(readability-braces-around-statements)
  if (spdlog::default_logger_raw()->name() == LOGGER_NAME) return;

  auto sink = std::make_shared<Sink>(std::forward<Args>(args)...);
  if constexpr (std::is_same_v<spdlog::async_logger, Logger>) {
    constexpr usize qsize = 32768;
    spdlog::init_thread_pool(qsize, 1);
    constexpr auto policy = spdlog::async_overflow_policy::block;
    spdlog::set_default_logger(std::make_shared<Logger>(LOGGER_NAME, std::move(sink), spdlog::thread_pool(), policy));
  } else {
    spdlog::set_default_logger(std::make_shared<Logger>(LOGGER_NAME, std::move(sink)));
  }

  spdlog::default_logger_raw()->flush_on(spdlog::level::err);
  spdlog::default_logger_raw()->set_pattern("%^%Y-%m-%dT%H:%M:%S.%e | [%L] | %v%$");
  spdlog::default_logger_raw()->set_level(spdlog::level::info);
}

inline void SetLancetLoggerLevel(spdlog::level::level_enum log_level) {
  spdlog::default_logger_raw()->set_level(log_level);
}

// NOLINTBEGIN(cppcoreguidelines-macro-usage)
#ifdef LANCET_DEVELOP_MODE
#define LOG_TRACE(...) spdlog::default_logger_raw()->trace(__VA_ARGS__);
#else
#define LOG_TRACE(...) ((void)0);
#endif
#define LOG_DEBUG(...) spdlog::default_logger_raw()->debug(__VA_ARGS__);
#define LOG_INFO(...) spdlog::default_logger_raw()->info(__VA_ARGS__);
#define LOG_WARN(...) spdlog::default_logger_raw()->warn(__VA_ARGS__);
#define LOG_ERROR(...) spdlog::default_logger_raw()->error(__VA_ARGS__);
#define LOG_CRITICAL(...) spdlog::default_logger_raw()->critical(__VA_ARGS__);
// NOLINTEND(cppcoreguidelines-macro-usage)

}  // namespace lancet

#endif  // SRC_LANCET_BASE_LOGGING_H_
