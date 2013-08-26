
#include <cassert>

#include "matrix.hh"

matrix::matrix(std::size_t rows, std::size_t columns) : rows_(rows), columns_(columns), data_(rows*columns) {
	// cout << "initializing a matrix..." << endl;
}
void matrix::set(std::size_t r, std::size_t c, double val) {
	data_[columns_*r + c] = val;
}
double matrix::get(std::size_t r, std::size_t c) const {
	return data_[columns_*r + c];
}
// double
std::ostream &matrix::print(std::ostream& os) const {
	for (std::size_t r=0; r<rows_; r++) {
		for (std::size_t c=0; c<columns_; c++) {
			if (c!=0) os << ", ";
			os << data_[columns_*r + c];
		}
		if (r+1 != rows_) os << "; ";
	}
    return os;
}
std::vector<double> matrix::multiply(const std::vector<double> &x) const {
	assert(x.size() == columns_);
	std::vector<double> retval(rows_);
	for (std::size_t current_row=0; current_row<rows_; ++current_row) {
		for (std::size_t current_col=0; current_col<columns_; ++current_col) {
			retval[current_row] += data_[columns_*current_row + current_col]*x[current_col];
		}
	}
	return retval;
}
std::vector<double> matrix::multiply_transpose(const std::vector<double> &y) const {
	assert(y.size() == rows_);
	// if (y.size() != rows_) throw 20;
	std::vector<double> retval(columns_);
	for (std::size_t current_row=0; current_row<rows_; ++current_row) {
		for (std::size_t current_col=0; current_col<columns_; ++current_col) {
			retval[current_col] += data_[columns_*current_row + current_col]*y[current_row];
		}
	}
	return retval;
}

