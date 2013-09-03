
#include <iostream>

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

sparse_matrix::sparse_matrix(std::size_t num_rows, std::size_t num_cols, std::size_t num_nonzeros) 
            : rows_(num_rows)
            , columns_(num_cols)
            , vals_(num_nonzeros)
            , row_indices_(num_nonzeros)
            , col_starts_(num_cols+1) {
    
}

sparse_matrix::sparse_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign) 
            : rows_(relevant_variables.size())
            , columns_(num_cols)
            , vals_(relevant_variables.size())
            , row_indices_(relevant_variables.size())
            , col_starts_(num_cols+1) {
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


sparse_matrix::sparse_matrix(std::size_t num_variables, double scalar) : rows_(num_variables), columns_(num_variables), vals_(num_variables), row_indices_(num_variables), col_starts_(num_variables+1) {
    if (num_variables == 0) return;
    
    for (size_t current_nonzero=0; current_nonzero < num_variables; ++current_nonzero) {
        vals_[current_nonzero] = scalar;
        row_indices_[current_nonzero] = current_nonzero;
    }
    for (size_t current_column=0; current_column < num_variables + 1; ++current_column) {
        col_starts_[current_column] = current_column;
    }
}

std::vector<double> sparse_matrix::multiply(const std::vector<double> &vector_factor) {
    assert(vector_factor.size() == columns_);
    std::vector<double> product(rows_);
    int current_row=-1;
    
    for (size_t column_index=0; column_index < columns_; ++column_index) {
        for (size_t nonzero_index=col_starts_[column_index]; nonzero_index < col_starts_[column_index+1]; ++nonzero_index) {
            current_row = row_indices_[nonzero_index];
            product[current_row] += vals_[nonzero_index]*vector_factor[column_index];
        }
    }
    return product;
}
std::vector<double> sparse_matrix::multiply_transpose(const std::vector<double> &vector_factor) {
    assert(vector_factor.size() == rows_);
    std::vector<double> product(columns_);
    int current_row=-1;
    
    for (size_t column_index=0; column_index < columns_; ++column_index) {
        for (size_t nonzero_index=col_starts_[column_index]; nonzero_index < col_starts_[column_index+1]; ++nonzero_index) {
            current_row = row_indices_[nonzero_index];
            product[column_index] += vals_[nonzero_index]*vector_factor[current_row];
        }
    }
    return product;
}

sparse_matrix vertical(const sparse_matrix &top, const sparse_matrix &bottom) {
    assert(top.num_columns() == bottom.num_columns());
    sparse_matrix result(   top.num_rows() + bottom.num_rows(), 
                            top.num_columns(), 
                            top.num_nnz() + bottom.num_nnz());
    size_t current_a_index=0, current_b_index=0;

    size_t current_nnz_index=0;
    for (size_t col_start_index=0; col_start_index < top.num_columns(); ++col_start_index) {
        // std::cout << "col_start_index: " <<col_start_index << std::endl;
        // std::cout << "              a: " << lower_var.col_starts_[col_start_index] << " -> " << lower_var.col_starts_[col_start_index+1] << std::endl;
        // std::cout << "              b: " << upper_var.col_starts_[col_start_index] << " -> " << upper_var.col_starts_[col_start_index+1] << std::endl;
        for (std::size_t current_a_index=top.col_starts_[col_start_index]; current_a_index < top.col_starts_[col_start_index+1]; ++current_a_index ) {
            // std::cout << "reading an A" << std::endl;
            result.vals_[current_nnz_index] = top.vals_[current_a_index];
            result.row_indices_[current_nnz_index] = top.row_indices_[current_a_index];
            ++current_nnz_index;
        }
        for (std::size_t current_b_index=bottom.col_starts_[col_start_index]; current_b_index < bottom.col_starts_[col_start_index+1]; ++current_b_index ) {
            // std::cout << "reading a B..." << std::endl;
            result.vals_[current_nnz_index] = bottom.vals_[current_b_index];
            result.row_indices_[current_nnz_index] = top.num_rows() + bottom.row_indices_[current_b_index];
            ++current_nnz_index;
        }
        result.col_starts_[col_start_index+1] = current_nnz_index;
    }

    // std::cout << "result:"<< std::endl << result << std::endl;
    return result;
}

sparse_matrix horizontal(const sparse_matrix &left, const sparse_matrix &right) {
    assert(left.num_rows() == right.num_rows());
    sparse_matrix result(   left.num_rows(),
                            left.num_columns() + right.num_columns(), 
                            left.num_nnz() + right.num_nnz());
    size_t result_nonzero_index=0;
    for (size_t left_nonzeros=0; left_nonzeros < left.num_nnz(); ++left_nonzeros) {
        result.vals_[result_nonzero_index] = left.vals_[left_nonzeros];
        result.row_indices_[result_nonzero_index] = left.row_indices_[left_nonzeros];
        ++result_nonzero_index;
    }
    for (size_t right_nonzeros=0; right_nonzeros < right.num_nnz(); ++right_nonzeros) {
        result.vals_[result_nonzero_index] = right.vals_[right_nonzeros];
        result.row_indices_[result_nonzero_index] = right.row_indices_[right_nonzeros];
        ++result_nonzero_index;
    }
    
    size_t result_column_index=0;
    for (size_t left_cols=0; left_cols < left.num_columns()+1; ++left_cols) {
        result.col_starts_[result_column_index] = left.col_starts_[left_cols];
        ++result_column_index;
    }
    result_column_index=left.num_columns();
    for (size_t right_cols=0; right_cols < right.num_columns()+1; ++right_cols) {
        result.col_starts_[result_column_index] = left.num_nnz() + right.col_starts_[right_cols];
        ++result_column_index;
    }
    return result;
}

