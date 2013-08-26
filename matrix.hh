
#ifndef __GUARD_MATRIX_HH
#define __GUARD_MATRIX_HH

#include <iostream>
#include <vector>


class matrix {
public:
	matrix(std::size_t rows, std::size_t columns);
	void set(std::size_t r, std::size_t c, double val);
	double get(std::size_t r, std::size_t c) const;
	// double
	std::ostream &print(std::ostream& os) const;
	std::vector<double> multiply(const std::vector<double> &x) const;
	std::vector<double> multiply_transpose(const std::vector<double> &y) const;
	std::size_t rows_, columns_;
	std::vector<double> data_;
private:
protected:
	
};

inline std::ostream& operator<< (std::ostream& os, const matrix& m) {
    os << "[";
	m.print(os);
    os << " ]";
    return os;
}

#endif
