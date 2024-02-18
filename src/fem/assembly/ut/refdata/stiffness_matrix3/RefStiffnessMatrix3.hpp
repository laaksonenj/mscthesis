#pragma once

#include "fem/multiprecision/Types.hpp"

namespace fem::ut::refdata
{
inline MatrixXmpq refStiffnessMatrix3 = (
	[]()
	{
		MatrixXmpq m(5,5);
		m(0,0) = mpq_class("1/2"); // 0.5
		m(0,1) = mpq_class("-1/2"); // -0.5
		m(0,2) = mpq_class("0"); // 0
		m(0,3) = mpq_class("0"); // 0
		m(0,4) = mpq_class("0"); // 0
		m(1,0) = mpq_class("-1/2"); // -0.5
		m(1,1) = mpq_class("3/2"); // 1.5
		m(1,2) = mpq_class("0"); // 0
		m(1,3) = mpq_class("0"); // 0
		m(1,4) = mpq_class("-1"); // -1
		m(2,0) = mpq_class("0"); // 0
		m(2,1) = mpq_class("0"); // 0
		m(2,2) = mpq_class("1"); // 1
		m(2,3) = mpq_class("0"); // 0
		m(2,4) = mpq_class("-1"); // -1
		m(3,0) = mpq_class("0"); // 0
		m(3,1) = mpq_class("0"); // 0
		m(3,2) = mpq_class("0"); // 0
		m(3,3) = mpq_class("1/2"); // 0.5
		m(3,4) = mpq_class("-1/2"); // -0.5
		m(4,0) = mpq_class("0"); // 0
		m(4,1) = mpq_class("-1"); // -1
		m(4,2) = mpq_class("-1"); // -1
		m(4,3) = mpq_class("-1/2"); // -0.5
		m(4,4) = mpq_class("5/2"); // 2.5
		return m;
	}
)();
} // namespace fem::ut::refdata
