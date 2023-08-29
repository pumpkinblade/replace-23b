#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "log.h"

namespace replace
{
  std::shared_ptr<spdlog::logger> Log::slogger_;

  void Log::Init()
  {
    std::vector<spdlog::sink_ptr> logSinks;
    logSinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    logSinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>("replace.log", true));

    logSinks[0]->set_pattern("%^[%T] [%l]: %v%$");
    logSinks[1]->set_pattern("[%T] [%l]: %v");

    slogger_ = std::make_shared<spdlog::logger>("replace", begin(logSinks), end(logSinks));
    spdlog::register_logger(slogger_);
    slogger_->set_level(spdlog::level::info);
    slogger_->flush_on(spdlog::level::info);
  }

}
