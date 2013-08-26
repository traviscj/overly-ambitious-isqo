#include <iostream>
#include <string>
#include <vector>
#include <cmath>

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

