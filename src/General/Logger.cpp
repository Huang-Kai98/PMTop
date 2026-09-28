//
// Created by huangkai on 9/17/25.
//

#include "Logger.h"
#include <vector>
#include <iostream>

Logger& Logger::getInstance() {
    static Logger instance;
    return instance;
}

void Logger::init(const std::string& logFile,
                  spdlog::level::level_enum level,
                  bool console,
                  bool file,
                  bool overwrite,
                  bool rotating,
                  size_t maxFileSize,
                  size_t maxFiles) {
    try {
        std::vector<spdlog::sink_ptr> sinks;

        if (console) {
            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            sinks.push_back(console_sink);
        }

        if (file) {
            if (rotating) {
                auto rotating_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
                    logFile, maxFileSize, maxFiles);
                sinks.push_back(rotating_sink);
            } else {
                auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                    logFile, overwrite);
                sinks.push_back(file_sink);
            }
        }

        logger_ = std::make_shared<spdlog::logger>("multi_sink", sinks.begin(), sinks.end());
        logger_->set_level(level);
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

        spdlog::register_logger(logger_);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "Log initialization failed: " << ex.what() << std::endl;
    }
}
