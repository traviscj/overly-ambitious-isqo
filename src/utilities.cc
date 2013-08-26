#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

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

