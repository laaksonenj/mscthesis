#pragma once

#include <any>
#include <cstdlib>
#include <functional>
#include <map>
#include <string>

#include <boost/program_options.hpp>

namespace fem
{
class Arguments
{
public:
    explicit Arguments(int argc, char* argv[]);

    template<typename T>
    T getValue(std::string option)
    {
        std::any value = m_optionParsers.at(option)(m_vm);
        if (!value.has_value())
        {
            std::exit(1);
        }
        return std::any_cast<T>(value);
    }

private:
    using OptionParserFn = std::function<std::any(const boost::program_options::variables_map&)>;

private:
    boost::program_options::options_description m_desc;
    boost::program_options::variables_map m_vm;
    std::map<std::string, OptionParserFn> m_optionParsers;
};
} // namespace fem
