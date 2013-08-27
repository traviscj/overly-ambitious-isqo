
#ifndef __GUARD_UTILITIES_HH
#define __GUARD_UTILITIES_HH

#include <iostream>
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

class sparse_matrix {
public:
    sparse_matrix(std::size_t num_rows, std::size_t num_cols, std::size_t num_nonzeros) : rows_(num_rows), columns_(num_cols), vals_(num_nonzeros), row_indices_(num_nonzeros), col_starts_(num_cols+1) {
        
    }
    
    // construct a 'relevant identity' matrix...
    sparse_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign) : rows_(relevant_variables.size()), columns_(num_cols), vals_(relevant_variables.size()), row_indices_(relevant_variables.size()), col_starts_(num_cols+1) {
        // sparse_matrix retval(relevant_variables.size(),total_variables);
        // double cur_sign = -1.0;
        std::size_t total_variable_index=0;
        for (std::size_t relevant_variable_index=0; relevant_variable_index < relevant_variables.size(); ++relevant_variable_index) {
            while (total_variable_index < relevant_variables[relevant_variable_index]) {
        
                ++total_variable_index;
                col_starts_[total_variable_index] = relevant_variable_index;
            }
            vals_[relevant_variable_index] = cur_sign;
            row_indices_[relevant_variable_index] = relevant_variable_index;
        }
        for (std::size_t i=total_variable_index; i < num_cols; ++i) {
            col_starts_[i+1] = col_starts_[total_variable_index]+1;
        }
    }
    
    
    std::size_t num_columns() const { return columns_;}
    std::size_t num_rows() const { return rows_;}
    std::size_t num_nnz() const { return vals_.size();}
    std::size_t rows_;
    std::size_t columns_;
    std::vector<double> vals_;
    std::vector<int> row_indices_;
    std::vector<int> col_starts_;
};

inline std::ostream& operator<< (std::ostream& os, const sparse_matrix& mat) {
    os << "rows: " << mat.num_rows() << "; cols: " << mat.num_columns() << "; nnz: " << mat.num_nnz() << std::endl;
    os << "vals       : " << mat.vals_ << std::endl;
    os << "row_indices: " << mat.row_indices_ << std::endl;
    os << "col_starts : " << mat.col_starts_ << std::endl;
    return os;
}

sparse_matrix vertical(const sparse_matrix &lower_var, const sparse_matrix &upper_var);

#endif
