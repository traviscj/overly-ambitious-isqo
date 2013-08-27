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

sparse_matrix vertical(const sparse_matrix &lower_var, const sparse_matrix &upper_var) {
    sparse_matrix result(   lower_var.num_rows() + upper_var.num_rows(), 
                            lower_var.num_columns(), 
                            lower_var.num_nnz() + upper_var.num_nnz());
    size_t current_a_index=0, current_b_index=0;

    size_t current_nnz_index=0;
    for (size_t col_start_index=0; col_start_index < lower_var.num_columns(); ++col_start_index) {
        // std::cout << "col_start_index: " <<col_start_index << std::endl;
        // std::cout << "              a: " << lower_var.col_starts_[col_start_index] << " -> " << lower_var.col_starts_[col_start_index+1] << std::endl;
        // std::cout << "              b: " << upper_var.col_starts_[col_start_index] << " -> " << upper_var.col_starts_[col_start_index+1] << std::endl;
        for (std::size_t current_a_index=lower_var.col_starts_[col_start_index]; current_a_index < lower_var.col_starts_[col_start_index+1]; ++current_a_index ) {
            // std::cout << "reading an A" << std::endl;
            result.vals_[current_nnz_index] = lower_var.vals_[current_a_index];
            result.row_indices_[current_nnz_index] = lower_var.row_indices_[current_a_index];
            ++current_nnz_index;
        }
        for (std::size_t current_b_index=upper_var.col_starts_[col_start_index]; current_b_index < upper_var.col_starts_[col_start_index+1]; ++current_b_index ) {
            // std::cout << "reading a B..." << std::endl;
            result.vals_[current_nnz_index] = upper_var.vals_[current_b_index];
            result.row_indices_[current_nnz_index] = lower_var.num_rows() + upper_var.row_indices_[current_b_index];
            ++current_nnz_index;
        }
        result.col_starts_[col_start_index+1] = current_nnz_index;
    }

    // std::cout << "result:"<< std::endl << result << std::endl;
    return result;
}