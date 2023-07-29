#ifndef __REPLACE_LOG__
#define __REPLACE_LOG__

#include <memory>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/ostr.h>

namespace replace
{
  class Log
  {
  public:
    static void Init();

    static std::shared_ptr<spdlog::logger> &getLogger() { return slogger_; }

  private:
    static std::shared_ptr<spdlog::logger> slogger_;
  };

}

// log macros
#define LOG_TRACE(...) ::replace::Log::getLogger()->trace(__VA_ARGS__)
#define LOG_DEBUG(...) ::replace::Log::getLogger()->debug(__VA_ARGS__)
#define LOG_INFO(...) ::replace::Log::getLogger()->info(__VA_ARGS__)
#define LOG_WARN(...) ::replace::Log::getLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...) ::replace::Log::getLogger()->error(__VA_ARGS__)
#define LOG_CRITICAL(...) ::replace::Log::getLogger()->critical(__VA_ARGS__)

#endif