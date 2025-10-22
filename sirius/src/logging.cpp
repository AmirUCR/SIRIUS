#include <iostream>
#include <cstdlib>     // std::exit

#include "utils.hpp"
#include "logging.hpp"

namespace sirius
{
    Logger logg{};

    void replace_all(std::string& source, const std::string& from, const std::string& to) 
    {
        std::string::size_type pos = 0;
        while ((pos = source.find(from, pos)) != std::string::npos) {
            source.replace(pos, from.length(), to);
            pos += to.length(); // Advance past the replacement
        }
    }

    void Logger::print_info(const std::string& message) const 
    {
        if (!quiet) std::cout << BLUE << "> " << RESET << message;
    }
    void Logger::print_info_newline(const std::string& message) const 
    {
        if (!quiet) std::cout << BLUE << "> " << RESET << message << '\n';
    }
    void Logger::print_timed_info_newline(const std::string& message) const
    {
        if (!quiet) std::cout << BLUE << "> " << RESET << "[" << utils::elapsed_since_start() << "] " << message << '\n';
    }
    void Logger::print_warning(const std::string& message) const
    {
        if (!quiet) std::cout << ORANGE << "> " << RESET << message;
    }
    void Logger::print_warning_newline(const std::string& message) const 
    {
        if (!quiet) std::cout << ORANGE << "> " << RESET << message << '\n';
    }
    void Logger::print_timed_warning_newline(const std::string& message) const
    {
        if (!quiet) std::cout << ORANGE << "> " << RESET << "[" << utils::elapsed_since_start() << "] " << message << '\n';
    }
    void Logger::print_error(const std::string& message) const 
    {
        if (!quiet) std::cout << RED << "> " << RESET << message;
    }
    void Logger::print_error_newline(const std::string& message) const 
    {
        if (!quiet) std::cout << RED << "> " << RESET << message << '\n';
    }
    void Logger::print_timed_error_newline(const std::string& message) const
    {
        if (!quiet) std::cout << RED << "> " << RESET << "[" << utils::elapsed_since_start() << "] " << message << '\n';
    }
    void Logger::begruessung() const 
    {
        if (!quiet) print_info_newline(std::string("Fly to ") + BLUE + "SIRIUS" + RESET);
    }

    void Logger::print_summary(const std::string& message) 
    {
        this->summary = message;
        
        if (!quiet) {
            std::string replacement = std::string("\n") + BLUE + ">" + RESET + " ";
            std::string mutable_message = message;

            replace_all(mutable_message, "\n", replacement);
            
            if (!mutable_message.empty() && mutable_message[0] == '\n') {
                mutable_message.replace(0, 1, "");
            }
            
            mutable_message = replacement + mutable_message + "\n";
            
            std::cout << mutable_message << std::endl;
        }
    }

    std::string Logger::format_error(const std::string& message) const 
    {
        return std::string(RED) + "> " + RESET + message;
    }

    [[noreturn]] void Logger::throw_formatted_error(const std::string& message) const 
    {
        if (!quiet) {
            std::cerr << format_error(message) << '\n';
            std::cout << format_error("Back to Earth.") << '\n';
        }
        std::exit(1);
    }

    [[noreturn]] void Logger::on_sigint()
    {
        if (!quiet) std::cout << "\n";
        print_error_newline("Interrupting SIRIUS. Ciao.");
        std::exit(1);
    }
} // namespace sirius
