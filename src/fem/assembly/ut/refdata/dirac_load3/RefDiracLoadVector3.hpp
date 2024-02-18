#pragma once

#include "fem/multiprecision/Types.hpp"

namespace fem::ut::refdata
{
inline VectorXmpq refDiracLoadVector3 = (
	[]()
	{
		VectorXmpq v(18);
		v(0) = mpq_class("0"); // 0
		v(1) = mpq_class("0"); // 0
		v(2) = mpq_class("1/4"); // 0.25
		v(3) = mpq_class("1/4"); // 0.25
		v(4) = mpq_class("1/2"); // 0.5
		v(5) = mpq_class("-1/4"); // -0.25
		v(6) = mpq_class("0"); // 0
		v(7) = mpq_class("-1/2"); // -0.5
		v(8) = mpq_class("-3/8"); // -0.375
		v(9) = mpq_class("-1/2"); // -0.5
		v(10) = mpq_class("-3/8"); // -0.375
		v(11) = mpq_class("0"); // 0
		v(12) = mpq_class("0"); // 0
		v(13) = mpq_class("0"); // 0
		v(14) = mpq_class("0"); // 0
		v(15) = mpq_class("0"); // 0
		v(16) = mpq_class("0"); // 0
		v(17) = mpq_class("1/32"); // 0.03125
		return v;
	}
)();
} // namespace fem::ut::refdata
