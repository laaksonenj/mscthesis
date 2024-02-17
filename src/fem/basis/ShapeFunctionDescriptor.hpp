#pragma once

#include <cstdint>
#include <tuple>
#include <variant>

#include <boost/functional/hash.hpp>

namespace fem
{
struct NodalShapeFunctionDescriptor
{
    uint32_t nodeIdx;
};

struct SideShapeFunctionDescriptor
{
    uint32_t sideIdx;
    uint32_t k;
};

struct InternalShapeFunctionDescriptor
{
    uint32_t k;
    uint32_t l;
};

using ShapeFunctionDescriptor = std::variant<NodalShapeFunctionDescriptor, SideShapeFunctionDescriptor, InternalShapeFunctionDescriptor>;

inline bool operator==(const NodalShapeFunctionDescriptor& lhs, const NodalShapeFunctionDescriptor& rhs)
{
    return lhs.nodeIdx == rhs.nodeIdx;
}

inline bool operator==(const SideShapeFunctionDescriptor& lhs, const SideShapeFunctionDescriptor& rhs)
{
    return lhs.sideIdx == rhs.sideIdx && lhs.k == rhs.k;
}

inline bool operator==(const InternalShapeFunctionDescriptor& lhs, const InternalShapeFunctionDescriptor& rhs)
{
    return lhs.k == rhs.k && lhs.l == rhs.l;
}

inline size_t hash_value(const NodalShapeFunctionDescriptor& val)
{
    return boost::hash_value(val.nodeIdx);
}

inline size_t hash_value(const SideShapeFunctionDescriptor& val)
{
    return boost::hash_value(std::make_pair(val.sideIdx, val.k));
}

inline size_t hash_value(const InternalShapeFunctionDescriptor& val)
{
    return boost::hash_value(std::make_pair(val.k, val.l));
}
} // namespace fem

namespace std
{
template<> struct hash<fem::NodalShapeFunctionDescriptor> : boost::hash<fem::NodalShapeFunctionDescriptor> {};
template<> struct hash<fem::SideShapeFunctionDescriptor> : boost::hash<fem::SideShapeFunctionDescriptor> {};
template<> struct hash<fem::InternalShapeFunctionDescriptor> : boost::hash<fem::InternalShapeFunctionDescriptor> {};
} // namespace std
