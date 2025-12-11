#ifndef _MULTI_THREAD_PRINTER_HPP
#define _MULTI_THREAD_PRINTER_HPP

/// 监控器，分文件记录组日志
/// @借助log4cp
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/consoleappender.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/config.hxx>
#include <log4cplus/configurator.h>
#include <log4cplus/spi/loggingevent.h>
#include <memory>
#include <iomanip>

// 单例
#define MultiThreadPrinter (MultiThreadPrinter::GetInstance())

// 颜色
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"

// 是否打印到控制台
#define ON_TERMINAL
// else ON_FIFE

inline char separator();

const char *file_name(const char *path);

#define REGISTER(key) (MultiThreadPrinter.Register(key))
#define INFO(key, msg) (MultiThreadPrinter.Info(key, msg, file_name(__FILE__), __LINE__, __func__))
#define DEBUG(key, msg) (MultiThreadPrinter.Debug(key, msg, file_name(__FILE__), __LINE__, __func__))
#define WARN(key, msg) (MultiThreadPrinter.Warn(key, msg, file_name(__FILE__), __LINE__, __func__))
#define ERROR(key, msg) (MultiThreadPrinter.Error(key, msg, file_name(__FILE__), __LINE__, __func__))

class MultiThreadPrinter {
private:
  std::map<std::string, log4cplus::Logger> monitors_t{};
  std::map<std::string, log4cplus::Logger> monitors_f{};

public:
  static MultiThreadPrinter &GetInstance() {
    static MultiThreadPrinter instance;
    return instance;
  }

  /// 对象构造时，将代号注册到日志记录器
  void Register(const std::string &key) {
    if (monitors_t.contains(key) || monitors_f.contains(key)) {
      return;
    }

#ifdef ON_TERMINAL
    monitors_t[key] = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT(key));
    const auto appender_t = log4cplus::SharedAppenderPtr(new log4cplus::ConsoleAppender());
    appender_t->setLayout(
      std::make_unique<log4cplus::PatternLayout>(LOG4CPLUS_TEXT("%d{%m/%d/%y %H:%M:%S.%q} [%t] %c - %m%n")));
    monitors_t[key].addAppender(appender_t);
    monitors_t[key].setAdditivity(false); // 禁止传播到根记录器
#else
        monitors_f[key] = log4cplus::Logger::getInstance(LOG4CPLUS_TEXT(key));
        const auto appender_f = log4cplus::SharedAppenderPtr(new log4cplus::FileAppender(LOG4CPLUS_TEXT("./log/" + key + ".log")));
        appender_f->setLayout(std::make_unique<log4cplus::PatternLayout>(LOG4CPLUS_TEXT("%d{%m/%d/%y %H:%M:%S.%q} [%t] %c - %m%n")));
        monitors_f[key].addAppender(appender_f);
        monitors_f[key].setAdditivity(false);  // 禁止日志传播到根记录器
#endif
  }

  void Info(const std::string &key, const std::string &msg,
            const char *file, int line, const char *func) {
    std::string location = " (" + std::string(file) + ":" + std::to_string(line) + " " + func + ")";
    assert(monitors_t.contains(key) || monitors_f.contains(key));
#ifdef ON_TERMINAL
    LOG4CPLUS_DEBUG(monitors_t[key], GREEN << msg << RESET << location);
#else
        LOG4CPLUS_DEBUG(monitors_f[key], msg << location);
#endif
  }

  void Debug(const std::string &key, const std::string &msg,
             const char *file, int line, const char *func) {
    std::string location = " (" + std::string(file) + ":" + std::to_string(line) + " " + func + ")";
    assert(monitors_t.contains(key) || monitors_f.contains(key));
#ifdef ON_TERMINAL
    LOG4CPLUS_DEBUG(monitors_t[key], CYAN << msg << RESET << location);
#else
        LOG4CPLUS_DEBUG(monitors_f[key], msg << location);
#endif
  }

  void Warn(const std::string &key, const std::string &msg,
            const char *file, int line, const char *func) {
    std::string location = " (" + std::string(file) + ":" + std::to_string(line) + " " + func + ")";
    assert(monitors_t.contains(key) || monitors_f.contains(key));
#ifdef ON_TERMINAL
    LOG4CPLUS_WARN(monitors_t[key], YELLOW << msg << RESET << location);
#else
        LOG4CPLUS_WARN(monitors_f[key], msg << location);
#endif
  }

  void Error(const std::string &key, const std::string &msg,
             const char *file, int line, const char *func) {
    std::string location = " (" + std::string(file) + ":" + std::to_string(line) + " " + func + ")";
    assert(monitors_t.contains(key) || monitors_f.contains(key));
#ifdef ON_TERMINAL
    LOG4CPLUS_ERROR(monitors_t[key], RED << msg << RESET << location);
#else
        LOG4CPLUS_ERROR(monitors_f[key], msg << location);
#endif
  }

  ~MultiThreadPrinter() {
    for (auto &librarian: monitors_t) {
      librarian.second.removeAllAppenders();
    }
    for (auto &librarian: monitors_f) {
      librarian.second.removeAllAppenders();
    }

    // 关闭日志记录器
    log4cplus::threadCleanup();
    log4cplus::Logger::shutdown();
  }

private:
  /// 构造
  MultiThreadPrinter() {
    // 初始化日志记录器
    log4cplus::initialize();
#ifndef ON_TERMINAL
        if (system("test -d ./log") != 0)
            system("mkdir -p ./log");
#endif
  }
};

inline char separator() {
#ifdef _WIN32
    return '\\';
#else
  return '/';
#endif
}

inline const char *file_name(const char *path) {
  const char *file = path;
  while (*path) {
    if (*path++ == separator()) {
      file = path;
    }
  }
  return file;
}

#endif // _MULTI_THREAD_PRINTER_HPP

