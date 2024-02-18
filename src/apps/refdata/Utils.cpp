#include "apps/refdata/Utils.hpp"

#include <sstream>

namespace fem
{
std::string generateCommonHeader()
{
    std::stringstream ss;
    ss << "#pragma once\n\n";
    ss << "#include \"fem/multiprecision/Types.hpp\"\n";
    return ss.str();
}

std::string namespaceBegin()
{
    std::stringstream ss;
    ss << "\nnamespace fem::ut::refdata\n";
    ss << "{\n";
    return ss.str();
}

std::string namespaceEnd()
{
    std::stringstream ss;
    ss << "} // namespace fem::ut::refdata\n";
    return ss.str();
}

std::string generateVariableDefinition(const MatrixXmpq& matrix, const std::string& variableName)
{
    std::stringstream ss;
    ss << "inline MatrixXmpq " << variableName << " = (\n";
    ss << "\t[]()\n";
    ss << "\t{\n";
    const uint32_t rows = matrix.rows();
    const uint32_t cols = matrix.cols();
    ss << "\t\tMatrixXmpq m(" << rows << "," << cols << ");\n";
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            const mpq_class v = matrix(row,col);
            ss << "\t\tm(" << row << "," << col << ") = mpq_class(\"" << v << "\")" << "; ";
            ss << "// " << v.get_d() << "\n";
        }
    }
    ss << "\t\treturn m;\n";
    ss << "\t}\n";
    ss << ")();\n";
    return ss.str();
}

std::string generateVariableDefinition(const VectorXmpq& vector, const std::string& variableName)
{
    std::stringstream ss;
    ss << "inline VectorXmpq " << variableName << " = (\n";
    ss << "\t[]()\n";
    ss << "\t{\n";
    const uint32_t dim = vector.size();
    ss << "\t\tVectorXmpq v(" << dim << ");\n";
    for (int i = 0; i < dim; i++)
    {
        const mpq_class v = vector(i);
        ss << "\t\tv(" << i << ") = mpq_class(\"" << v << "\")" << "; ";
        ss << "// " << v.get_d() << "\n";
    }
    ss << "\t\treturn v;\n";
    ss << "\t}\n";
    ss << ")();\n";
    return ss.str();
}
} // namespace fem
