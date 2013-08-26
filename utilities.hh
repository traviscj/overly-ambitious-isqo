
#ifndef __GUARD_UTILITIES_HH
#define __GUARD_UTILITIES_HH
#include <string>
#include <vector>

double bracket_plus(double val);
double bracket_minus(double val);
std::string ordinal(int n);
bool assert_close(double val1, double val2, double tol);

// want to print vector<int> and vector<double>...
template < class T >
inline std::ostream& operator<< (std::ostream& os, const std::vector< T >& vec) {
    os << "[";
	for (size_t vector_index = 0; vector_index < vec.size(); ++vector_index) {
		if (vector_index != 0) os << ", ";
		os << vec[vector_index];
	}
    os << " ]";
    return os;
}

#endif
