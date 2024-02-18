#pragma once

#include "fem/multiprecision/Types.hpp"

namespace fem::ut::refdata
{
inline VectorXmpq refDiracLoadVector1 = (
	[]()
	{
		VectorXmpq v(27);
		v(0) = mpq_class("0"); // 0
		v(1) = mpq_class("0"); // 0
		v(2) = mpq_class("0"); // 0
		v(3) = mpq_class("0"); // 0
		v(4) = mpq_class("1"); // 1
		v(5) = mpq_class("0"); // 0
		v(6) = mpq_class("0"); // 0
		v(7) = mpq_class("0"); // 0
		v(8) = mpq_class("0"); // 0
		v(9) = mpq_class("0"); // 0
		v(10) = mpq_class("0"); // 0
		v(11) = mpq_class("0"); // 0
		v(12) = mpq_class("0"); // 0
		v(13) = mpq_class("0"); // 0
		v(14) = mpq_class("0"); // 0
		v(15) = mpq_class("0"); // 0
		v(16) = mpq_class("0"); // 0
		v(17) = mpq_class("0"); // 0
		v(18) = mpq_class("0"); // 0
		v(19) = mpq_class("0"); // 0
		v(20) = mpq_class("0"); // 0
		v(21) = mpq_class("0"); // 0
		v(22) = mpq_class("0"); // 0
		v(23) = mpq_class("0"); // 0
		v(24) = mpq_class("0"); // 0
		v(25) = mpq_class("0"); // 0
		v(26) = mpq_class("0"); // 0
		return v;
	}
)();
} // namespace fem::ut::refdata
