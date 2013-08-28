
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
        if (relevant_variables.size()==0) return;
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
    
    // pure identity:
    // values_{I_n}: scalar*[1.0, 1.0, 1.0,...] (num_variables entries)
    // row_indices_{I_n}: [0,1,2,...,num_variables-1] (num_variables entries)
    // col_starts_{I_n}: [0,1,2,...,num_variables] (num_variables+1 entries)
    sparse_matrix(std::size_t num_variables, double scalar) : rows_(num_variables), columns_(num_variables), vals_(num_variables), row_indices_(num_variables), col_starts_(num_variables+1) {
        if (num_variables == 0) return;
        
        for (size_t current_nonzero=0; current_nonzero < num_variables; ++current_nonzero) {
            vals_[current_nonzero] = scalar;
            row_indices_[current_nonzero] = current_nonzero;
        }
        for (size_t current_column=0; current_column < num_variables + 1; ++current_column) {
            col_starts_[current_column] = current_column;
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

// TODO: 'add_to_diag' function for sparse_matrix.
// TODO: change this to 'append_below' & put into sparse_matrix class.
sparse_matrix vertical(const sparse_matrix &lower_var, const sparse_matrix &upper_var);
// TODO: change this to 'append_right' & put into sparse_matrix class.
sparse_matrix horizontal(const sparse_matrix &left, const sparse_matrix &right);

// TODO: move sparse_matrix class into it's own class.

#endif
