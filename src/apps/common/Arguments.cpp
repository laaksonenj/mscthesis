#include "apps/common/Arguments.hpp"

#include <iostream>

#include "fem/basis/PolynomialSpaceType.hpp"

namespace fem
{
namespace po = boost::program_options;

Arguments::Arguments(int argc, char* argv[])
{
    m_desc.add_options()
        ("mesh-file", po::value<std::string>(), "")
        ("p", po::value<int>(), "")
        ("polynomial-space", po::value<std::string>())
        ("variable-name", po::value<std::string>())
        ;

    po::store(po::parse_command_line(argc, argv, m_desc), m_vm);
    po::notify(m_vm);

    m_optionParsers.emplace("mesh-file", [](const po::variables_map& vm)
    {
        if (vm.count("mesh-file"))
        {
            return std::any(vm["mesh-file"].as<std::string>());
        }
        else
        {
            std::cout << "Argument missing: --mesh-file" << std::endl;
            return std::any();
        }
    });

    m_optionParsers.emplace("p", [](const po::variables_map& vm)
    {
        if (vm.count("p"))
        {
            const int p = vm["p"].as<int>();
            if (p >= 1)
            {
                return std::any(p);
            }
            else
            {
                std::cout << "p must be greater than or equal to 1" << std::endl;
                return std::any();
            }
        }
        else
        {
            std::cout << "Argument missing: --p" << std::endl;
            return std::any();
        }
    });

    m_optionParsers.emplace("polynomial-space", [](const po::variables_map& vm)
    {
        if (vm.count("polynomial-space"))
        {
            const std::string str = vm["polynomial-space"].as<std::string>();
            if (str == "product")
            {
                return std::any(PolynomialSpaceType_Product);
            }
            else if (str == "trunk")
            {
                return std::any(PolynomialSpaceType_Trunk);
            }
            else
            {
                std::cout << "Invalid polynomial space type: " << str << std::endl;
                return std::any();
            }
        }
        else
        {
            std::cout << "Argument missing: --polynomial-space" << std::endl;
            return std::any();
        }
    });

    m_optionParsers.emplace("variable-name", [](const po::variables_map& vm)
    {
        if (vm.count("variable-name"))
        {
            return std::any(vm["variable-name"].as<std::string>());
        }
        else
        {
            std::cout << "Argument missing: --variable-name" << std::endl;
            return std::any();
        }
    });
}
} // namespace fem
