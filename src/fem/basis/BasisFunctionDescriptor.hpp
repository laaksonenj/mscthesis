#pragma once

#include <cstdint>
#include <tuple>
#include <variant>

#include <boost/functional/hash.hpp>

#include "fem/domain/Mesh.hpp"

namespace fem
{
struct NodalBasisFunctionDescriptor
{
    Mesh::NodeIndex nodeIdx;
};

struct SideBasisFunctionDescriptor
{
    Mesh::SideIndex sideIdx;
    uint32_t k;
};

struct InternalBasisFunctionDescriptor
{
    Mesh::ElementIndex elementIdx;
    uint32_t k;
    uint32_t l;
};

using BasisFunctionDescriptor = std::variant<NodalBasisFunctionDescriptor, SideBasisFunctionDescriptor, InternalBasisFunctionDescriptor>;

inline bool operator==(const NodalBasisFunctionDescriptor& lhs, const NodalBasisFunctionDescriptor& rhs)
{
    return lhs.nodeIdx == rhs.nodeIdx;
}

inline bool operator==(const SideBasisFunctionDescriptor& lhs, const SideBasisFunctionDescriptor& rhs)
{
    return lhs.sideIdx == rhs.sideIdx && lhs.k == rhs.k;
}

inline bool operator==(const InternalBasisFunctionDescriptor& lhs, const InternalBasisFunctionDescriptor& rhs)
{
    return lhs.elementIdx == rhs.elementIdx && lhs.k == rhs.k && lhs.l == rhs.l;
}

inline size_t hash_value(const NodalBasisFunctionDescriptor& val)
{
    return boost::hash_value(val.nodeIdx);
}

inline size_t hash_value(const SideBasisFunctionDescriptor& val)
{
    return boost::hash_value(std::make_pair(val.sideIdx, val.k));
}

inline size_t hash_value(const InternalBasisFunctionDescriptor& val)
{
    return boost::hash_value(std::make_tuple(val.elementIdx, val.k, val.l));
}
} // namespace fem

namespace std
{
template<> struct hash<fem::NodalBasisFunctionDescriptor> : boost::hash<fem::NodalBasisFunctionDescriptor> {};
template<> struct hash<fem::SideBasisFunctionDescriptor> : boost::hash<fem::SideBasisFunctionDescriptor> {};
template<> struct hash<fem::InternalBasisFunctionDescriptor> : boost::hash<fem::InternalBasisFunctionDescriptor> {};
} // namespace std
