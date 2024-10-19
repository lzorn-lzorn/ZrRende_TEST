#pragma once

#include <mutex>
#include <string>
#include <string_view>
#include <vector>
#include <format>
#include <fstream>
#include <iostream>

namespace ZrRender {
#if defined(_MSVC_VER)
    #define __CALLER__    __FUNCSIG__
    #define __FILE_NAME__ __FILE_NAME__
#elif defined(__clang__) || defined(__GNUC__)
    #define __CALLER__ __PRETTY_FUNCTION__

#else
    #error "Unsupported Compiler"
#endif

std::string GetFileName(const std::string &);

#define NOTICE(fmt, ...) std::cout << std::format("[{}:{}] " fmt, __FILE_NAME__, __LINE__, __VA_ARGS__) << std::endl;

enum class LogLevel;
enum class LogTarget;
class Log;

namespace detail {

struct FileManager {};

struct Message {
    LogLevel    m_level;
    LogTarget   m_target;
    size_t      m_line;
    std::string m_caller;
    std::string m_file;
    std::string m_buffer;
};

}  // namespace detail
enum class LogLevel {
    NONE    = 0,
    INFO    = 1,
    DEBUG   = 2,
    WARNING = 3,
    ERROR   = 4,
};
enum class LogTarget {
    NONE    = 0,
    CONSOLE = 1,
    FILE    = 2,
    BOTH    = 3,
};

class Log {
private:
    const std::mutex             _p_mutex;
    std::vector<detail::Message> _p_messages;
    // detail::FileManager          _p_file;
    std::string       _p_file_path = "C:/Program User/code/ZrRender/LOG/";
    std::string       _p_file_name = "LOG.txt";
    const std::string _p_file_format;
    const std::string _p_console_format;

public:
    Log(const Log &)            = delete;
    Log(Log &&)                 = delete;
    Log &operator=(const Log &) = delete;
    Log &operator=(Log &&)      = delete;

    ~Log() {
        // std::lock_guard lock(_p_mutex);
        std::ofstream out_stream;
        out_stream.open(_p_file_path + _p_file_name, std::ios::app | std::ios::out);
        if (!out_stream.is_open()) {
            std::cerr << "Failed to open log file: " << _p_file_path + _p_file_name << std::endl;
            out_stream.close();
        } else {
            for (auto every_message : _p_messages) {
                auto log_context = every_message.m_buffer.c_str();
                out_stream.write(log_context, every_message.m_buffer.size());
            }
            out_stream.close();
        }
    }

    void SetLevel(LogLevel level) noexcept;
    void SetTarget(LogTarget target) noexcept;

    void SetFilePath(std::string path, std::string name) {
        _p_file_path = path + name;
    }

    template <typename... Args>
    void LOG(LogLevel level, LogTarget target, std::string format, Args &&...args) noexcept {
        std::string     head_format = "[{},{},level:{}]:Message:[\n\t your_contexts \n]\n ";
        detail::Message message;
        message.m_level      = level;
        message.m_target     = target;
        message.m_line       = __LINE__;
        message.m_caller     = __CALLER__;
        message.m_file       = __FILE_NAME__;
        std::string tail_str = std::string(std::vformat(format, std::make_format_args(std::forward<Args>(args)...)));
        message.m_buffer =
            std::string(std::format("[{},{},level:{}]:Message:[ ", __FILE_NAME__, __LINE__, EnumToString(level)))
            + tail_str + "]\n";
        std::cout << message.m_buffer << std::endl;
        _p_messages.push_back(message);
    }

public:
    static Log &Instance() noexcept {
        static Log instance;
        return instance;
    }

private:
    Log() {
    }

    static Log *logger;

    std::string EnumToString(LogLevel level) const noexcept {
        switch (level) {
            case LogLevel::INFO:
                return "Info";
            case LogLevel::WARNING:
                return "Warning";
            case LogLevel::ERROR:
                return "Error";
            default:
                return "Unknown";
        }
    }
};

}  // namespace ZrRender