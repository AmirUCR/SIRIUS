#pragma once
#include <string>

#define RESET "\033[0m"
#define RED "\033[31m"
#define BLUE "\033[34;1m"
#define ORANGE "\033[38;5;208m"

namespace sirius
{
    struct Logger
    {
        bool quiet = false;
        std::string summary;

        void begruessung() const;
        void print_info(const std::string &message) const;
        void print_info_newline(const std::string &message) const;
        void print_timed_info_newline(const std::string& message) const;
        void print_warning(const std::string &message) const;
        void print_warning_newline(const std::string &message) const;
        void print_timed_warning_newline(const std::string& message) const;
        void print_error(const std::string &message) const;
        void print_error_newline(const std::string &message) const;
        void print_timed_error_newline(const std::string& message) const;

        void print_summary(const std::string &message);
        std::string format_error(const std::string &message) const;
        [[noreturn]] void throw_formatted_error(const std::string &message) const;
        [[noreturn]] void on_sigint();
    };

    extern Logger logg; // defined in logging.cpp
}
