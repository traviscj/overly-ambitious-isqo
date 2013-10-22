// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

#include "utilities.hh"

double bracket_plus(double val) {
	return fmax(val, 0.0);
}
double bracket_minus(double val) {
	return fmax(-val, 0.0);
}
std::string ordinal(int n) {
	if (n%10==1)
		return std::string("st");
	if (n%10==2)
		return std::string("nd");
	if (n%10==3)
		return std::string("rd");
	else
		return std::string("th");
}

bool assert_close(double val1, double val2, double tol) {
	assert(abs(val1) <= abs(val2) + tol*abs(val1));
	assert(abs(val2) <= abs(val1) + tol*abs(val2));
	assert(((val1 > tol) && (val2 > tol)) || ((val1 < -tol) && (val2 < -tol)) || (val1 == 0 && val2 == 0));
	return (abs(val1) <= abs(val2) + tol*abs(val1)) && (abs(val2) <= abs(val1) + tol*abs(val2));
}

double dot_product(const std::vector<double> &vector_one, const std::vector<double> &vector_two) {
    assert(vector_one.size() == vector_two.size());
	double dot_product=0.0;
	for (size_t primal_index=0; primal_index < vector_one.size(); ++primal_index) {
		dot_product += vector_one[primal_index]*vector_two[primal_index];
	}
    return dot_product;
}