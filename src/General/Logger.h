//
// Created by huangkai on 9/17/25.
//
#pragma once
#include <memory>
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>

class Logger {
public:
    static Logger& getInstance();

    /**
     * 初始化日志
     * @param logFile 日志文件名
     * @param level 日志级别
     * @param console 是否输出到控制台
     * @param file 是否输出到文件
     * @param overwrite 文件是否覆盖（仅 basic_file_sink）
     * @param rotating 是否轮转日志
     * @param maxFileSize 最大文件大小（轮转日志用）
     * @param maxFiles 最大文件个数（轮转日志用）
     */
    void init(const std::string& logFile = "log.txt",
              spdlog::level::level_enum level = spdlog::level::info,
              bool console = true,
              bool file = true,
              bool overwrite = false,
              bool rotating = false,
              size_t maxFileSize = 1048576,
              size_t maxFiles = 3);

    // 日志接口
    template<typename... Args>
    void info(const char* fmt, const Args&... args) { if (logger_) logger_->info(fmt, args...); }
    template<typename... Args>
    void warn(const char* fmt, const Args&... args) { if (logger_) logger_->warn(fmt, args...); }
    template<typename... Args>
    void error(const char* fmt, const Args&... args) { if (logger_) logger_->error(fmt, args...); }
    template<typename... Args>
    void debug(const char* fmt, const Args&... args) { if (logger_) logger_->debug(fmt, args...); }
    template<typename... Args>
    void trace(const char* fmt, const Args&... args) { if (logger_) logger_->trace(fmt, args...); }

private:
    Logger() = default;
    ~Logger() = default;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    std::shared_ptr<spdlog::logger> logger_;
};

#define LOG_INFO(fmt, ...)  Logger::getInstance().info("[{}:{}] " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define LOG_WARN(fmt, ...)  Logger::getInstance().warn("[{}:{}] " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define LOG_ERROR(fmt, ...) Logger::getInstance().error("[{}:{}] " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define LOG_DEBUG(fmt, ...) Logger::getInstance().debug("[{}:{}] " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
#define LOG_TRACE(fmt, ...) Logger::getInstance().trace("[{}:{}] " fmt, __FILE__, __LINE__, ##__VA_ARGS__)
