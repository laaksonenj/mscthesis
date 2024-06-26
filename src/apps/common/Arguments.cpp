#include "apps/common/Arguments.hpp"

#include <algorithm>
#include <iostream>

#include "apps/common/LinearSolver.hpp"
#include "fem/basis/PolynomialSpaceType.hpp"
#include "fem/multiprecision/Types.hpp"

namespace fem
{
namespace po = boost::program_options;

#define ARGUMENT_MISSING(argName) \
    std::cout << "Argument missing: --" << argName << std::endl; \
    return std::any();

Arguments::Arguments(int argc, char* argv[])
{
    m_desc.add_options()
        ("mesh-file", po::value<std::string>())
        ("p", po::value<int>())
        ("polynomial-space", po::value<std::string>())
        ("variable-name", po::value<std::string>())
        ("dirac-point", po::value<std::vector<std::string>>()->multitoken())
        ("boundary-function-idx", po::value<int>())
        ("element-idx", po::value<int>())
        ("local-side-idx", po::value<int>())
        ("output-dir", po::value<std::string>()->default_value(""))
        ("precision", po::value<int>()->default_value(64))
        ("linear-solver", po::value<std::string>()->default_value("ldlt"))
        ;

    po::store(po::parse_command_line(argc, argv, m_desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), m_vm);
    po::notify(m_vm);

    m_optionParsers.emplace("mesh-file", [](const po::variables_map& vm)
    {
        if (vm.count("mesh-file"))
        {
            return std::any(vm["mesh-file"].as<std::string>());
        }
        else
        {
            ARGUMENT_MISSING("mesh-file");
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
            ARGUMENT_MISSING("p");
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
            ARGUMENT_MISSING("polynomial-space");
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
            ARGUMENT_MISSING("variable-name");
        }
    });

    m_optionParsers.emplace("dirac-point", [](const po::variables_map& vm)
    {
        if (vm.count("dirac-point"))
        {
            const std::vector<std::string> coords = vm["dirac-point"].as<std::vector<std::string>>();
            if (coords.size() != 2)
            {
                std::cout << "Dirac point must contain two coordinate values" << std::endl;
                return std::any();
            }
            const Vector2mpq x_0{mpq_class(coords[0]), mpq_class(coords[1])};
            return std::any(x_0);
        }
        else
        {
            ARGUMENT_MISSING("dirac-point");
        }
    });

    m_optionParsers.emplace("boundary-function-idx", [](const po::variables_map& vm)
    {
        if (vm.count("boundary-function-idx"))
        {
            const int boundaryFunctionIdx = vm["boundary-function-idx"].as<int>();
            if (boundaryFunctionIdx >= 0)
            {
                return std::any(boundaryFunctionIdx);
            }
            else
            {
                return std::any();
            }
        }
        else
        {
            ARGUMENT_MISSING("boundary-function-idx");
        }
    });

    m_optionParsers.emplace("element-idx", [](const po::variables_map& vm)
    {
        if (vm.count("element-idx"))
        {
            const int elementIdx = vm["element-idx"].as<int>();
            if (elementIdx >= 0)
            {
                return std::any(elementIdx);
            }
            else
            {
                return std::any();
            }
        }
        else
        {
            ARGUMENT_MISSING("element-idx");
        }
    });

    m_optionParsers.emplace("local-side-idx", [](const po::variables_map& vm)
    {
        if (vm.count("local-side-idx"))
        {
            const int localSideIdx = vm["local-side-idx"].as<int>();
            if (localSideIdx >= 0)
            {
                return std::any(localSideIdx);
            }
            else
            {
                return std::any();
            }
        }
        else
        {
            ARGUMENT_MISSING("local-side-idx");
        }
    });

    m_optionParsers.emplace("output-dir", [](const po::variables_map& vm)
    {
        if (vm.count("output-dir"))
        {
            return std::any(vm["output-dir"].as<std::string>());
        }
        else
        {
            ARGUMENT_MISSING("output-dir");
        }
    });

    m_optionParsers.emplace("precision", [](const po::variables_map& vm)
    {
        if (vm.count("precision"))
        {
            const int precision = vm["precision"].as<int>();
            if (precision >= 64)
            {
                return std::any(static_cast<uint32_t>(precision));
            }
            else
            {
                std::cout << "precision must be greater than or equal to 64" << std::endl;
                return std::any();
            }
        }
        else
        {
            ARGUMENT_MISSING("precision");
        }
    });

    m_optionParsers.emplace("linear-solver", [](const po::variables_map& vm)
    {
        if (vm.count("linear-solver"))
        {
            const std::string method = vm["linear-solver"].as<std::string>();
            const auto it = std::find_if(linearSolverMethodCliNames.begin(), linearSolverMethodCliNames.end(),
                [&method](const auto& elem)
                {
                    return elem.second == method;
                }
            );
            if (it != linearSolverMethodCliNames.end())
            {
                return std::any(it->first);
            }
            else
            {
                std::cout << "Invalid linear solver method: " << method << std::endl;
                return std::any();
            }
        }
        else
        {
            ARGUMENT_MISSING("linear-solver");
        }
    });
}
} // namespace fem
