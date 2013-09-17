
#include <iostream>

#include <cassert>

#include "matrix.hh"

matrix_base_class::matrix_base_class(std::size_t rows, std::size_t columns) 
                    : rows_(rows)
                    , columns_(columns)
                    , last_hessian_shift_(-3.14159) { 
}


matrix_base_class::~matrix_base_class() {
    // std::cout << "destroying matrix base class" << std::endl;
}

dense_matrix::dense_matrix(std::size_t rows, std::size_t columns) : matrix_base_class(rows, columns), data_(rows*columns) {
	// cout << "initializing a matrix..." << endl;
}
dense_matrix::dense_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign) : matrix_base_class(num_cols, relevant_variables.size()), data_(rows_*columns_) {
    for (size_t current_relevant_variable=0; current_relevant_variable < relevant_variables.size(); ++current_relevant_variable) {
        set(current_relevant_variable, relevant_variables[current_relevant_variable], cur_sign);
    }
}    
dense_matrix::dense_matrix(std::size_t num_variables, double scalar) : matrix_base_class(num_variables, num_variables), data_(rows_*columns_) {
    for (size_t current_variable=0; current_variable < num_variables; ++current_variable) {
        set(current_variable, current_variable, scalar);
    }
}

void dense_matrix::set(std::size_t r, std::size_t c, double val) {
	data_[columns_*r + c] = val;
}
double dense_matrix::get(std::size_t r, std::size_t c) const {
	return data_[columns_*r + c];
}

std::ostream &dense_matrix::print(std::ostream& os) const {
	for (std::size_t r=0; r<rows_; r++) {
		for (std::size_t c=0; c<columns_; c++) {
			if (c!=0) os << ", ";
			os << data_[columns_*r + c];
		}
		if (r+1 != rows_) os << "; ";
	}
    return os;
}
std::ostream &sparse_matrix::print(std::ostream& os) const {
    os << "values       : " << this->vals_ << std::endl;
    os << "row indices  : " << this->row_indices_ << std::endl;
    os << "column starts: " << this->col_starts_ << std::endl;
    
    for (std::size_t col_index=0; col_index < this->columns_; ++col_index) {
        os << "col: " << col_index << "\t";
        for (std::size_t nnz_index=this->col_starts_[col_index]; nnz_index < this->col_starts_[col_index+1]; ++nnz_index) {
            if (nnz_index != this->col_starts_[col_index]) os << "\t";
            os << "row: " << this->row_indices_[nnz_index] << "; val:" << this->vals_[nnz_index] << std::endl;
        }
    }
    return os;
}

std::vector<double> dense_matrix::multiply(const std::vector<double> &x) const {
	assert(x.size() == columns_);
	std::vector<double> retval(rows_);
	for (std::size_t current_row=0; current_row<rows_; ++current_row) {
		for (std::size_t current_col=0; current_col<columns_; ++current_col) {
			retval[current_row] += data_[columns_*current_row + current_col]*x[current_col];
		}
	}
	return retval;
}
std::vector<double> dense_matrix::multiply_transpose(const std::vector<double> &y) const {
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
            : matrix_base_class(num_rows, num_cols)
            , vals_(num_nonzeros)
            , row_indices_(num_nonzeros)
            , col_starts_(num_cols+1) {
    
}

sparse_matrix::sparse_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign) 
            : matrix_base_class(relevant_variables.size(), num_cols)
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


sparse_matrix::sparse_matrix(std::size_t num_variables, double scalar)
            : matrix_base_class(num_variables, num_variables)
            , vals_(num_variables)
            , row_indices_(num_variables)
            , col_starts_(num_variables+1) {
    if (num_variables == 0) return;
    
    for (size_t current_nonzero=0; current_nonzero < num_variables; ++current_nonzero) {
        vals_[current_nonzero] = scalar;
        row_indices_[current_nonzero] = current_nonzero;
    }
    for (size_t current_column=0; current_column < num_variables + 1; ++current_column) {
        col_starts_[current_column] = current_column;
    }
}

std::vector<double> sparse_matrix::multiply(const std::vector<double> &vector_factor) const {
    assert(vector_factor.size() == columns_);
    std::vector<double> product(rows_);
    // std::cout << "product, pre: " << product << std::endl;
    for (size_t row_index=0; row_index < rows_; ++row_index) product[row_index]=0.0;
    int current_row=-1;
    
    for (size_t column_index=0; column_index < columns_; ++column_index) {
        for (size_t nonzero_index=col_starts_[column_index]; nonzero_index < col_starts_[column_index+1]; ++nonzero_index) {
            current_row = row_indices_[nonzero_index];
            product[current_row] += vals_[nonzero_index]*vector_factor[column_index];
        }
    }
    // std::cout << "product, post: " << product << std::endl;
    return product;
}
std::vector<double> sparse_matrix::multiply_transpose(const std::vector<double> &vector_factor) const {
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

void dense_matrix::sum(const std::shared_ptr<dense_matrix> b) {
    assert(num_rows() == b->num_rows());
    assert(num_columns() == b->num_columns());
    std::shared_ptr<dense_matrix> result(new dense_matrix(num_rows(), num_columns()));
    
    for (size_t row_index=0; row_index < num_rows(); ++row_index) {
        for (size_t column_index=0; column_index < num_rows(); ++column_index) {
            set(row_index, column_index, get(row_index, column_index) + b->get(row_index, column_index));
        }
    }
    // return result;
}
void sparse_matrix::sum(const std::shared_ptr<sparse_matrix> b) {
    assert(num_rows() == b->num_rows());
    assert(num_columns() == b->num_columns());
    bool PRINT = false;
    // std::shared_ptr<sparse_matrix> result(new sparse_matrix(num_rows(), a->num_columns(), 0));
    
    if (PRINT) std::cout << std::endl << "summing: " << std::endl << *this << std::endl;
    if (PRINT) std::cout << "and : " << std::endl << *b << std::endl;
    std::vector<double> result_vals;
    std::vector<int> result_row_indices;
    std::vector<int> result_col_starts;
    size_t A_nnz_index=0;
    size_t B_nnz_index=0;
    size_t C_nnz_index=0, C_col_index=0;
    result_col_starts.push_back(C_nnz_index);
    for (; C_col_index < col_starts_.size() - 1; ) {        
        if (PRINT) std::cout << "A bounds: " << col_starts_[C_col_index] << " <= A_nnz_index < " << col_starts_[C_col_index+1] << "; "
                             << "B bounds: " << b->col_starts_[C_col_index] << " <= B_nnz_index < " << b->col_starts_[C_col_index+1]
                             << std::endl;
        
        A_nnz_index=col_starts_[C_col_index];
        B_nnz_index=b->col_starts_[C_col_index];
        while ( (A_nnz_index < col_starts_[C_col_index+1]) &&  (B_nnz_index < b->col_starts_[C_col_index+1]) ) {
            if (row_indices_[A_nnz_index] < b->row_indices_[B_nnz_index]) {
                result_vals.push_back(vals_[A_nnz_index]);
                result_row_indices.push_back(row_indices_[A_nnz_index]);
                ++A_nnz_index;
                ++C_nnz_index;
            } else if (row_indices_[A_nnz_index] > b->row_indices_[B_nnz_index]) {
                result_vals.push_back(b->vals_[B_nnz_index]);
                result_row_indices.push_back(b->row_indices_[B_nnz_index]);
                ++B_nnz_index;
                ++C_nnz_index;
            } else {
                assert((row_indices_[A_nnz_index] == b->row_indices_[B_nnz_index]));
                result_vals.push_back(vals_[A_nnz_index] + b->vals_[B_nnz_index]);
                result_row_indices.push_back(row_indices_[A_nnz_index]);
                ++A_nnz_index;
                ++B_nnz_index;
                ++C_nnz_index;
            }
        }
                
        while ( (A_nnz_index < col_starts_[C_col_index+1]) ) {
            result_vals.push_back(vals_[A_nnz_index]);
            result_row_indices.push_back(row_indices_[A_nnz_index]);
            ++A_nnz_index;
            ++C_nnz_index;
        }
        while ( (B_nnz_index < b->col_starts_[C_col_index+1]) ) {
            result_vals.push_back(b->vals_[B_nnz_index]);
            result_row_indices.push_back(b->row_indices_[B_nnz_index]);
            ++B_nnz_index;
            ++C_nnz_index;
        }
        
        result_col_starts.push_back(C_nnz_index);
        ++C_col_index;
        
    }
    vals_ = result_vals;
    row_indices_ = result_row_indices;
    col_starts_ = result_col_starts;
    if (PRINT) std::cout << "result : " << std::endl << *this << std::endl << "=========" << std::endl;
    
}
// TODO move this into dense_matrix
std::shared_ptr<dense_matrix> vertical(const std::shared_ptr<dense_matrix> top, const std::shared_ptr<dense_matrix> bottom) {
    assert(top->num_columns() == bottom->num_columns());
    std::shared_ptr<dense_matrix> result(new dense_matrix(top->num_rows() + bottom->num_rows(), 
                            top->num_columns()));

    for (size_t current_top_row_index=0; current_top_row_index < top->num_rows(); ++current_top_row_index) {
        for (size_t current_top_col_index=0; current_top_col_index < top->num_columns(); ++current_top_col_index) {
            result->set(current_top_row_index, current_top_col_index, top->get(current_top_row_index, current_top_col_index));
        }
    }
    for (size_t current_bottom_row_index=0; current_bottom_row_index < bottom->num_rows(); ++current_bottom_row_index) {
        for (size_t current_bottom_col_index=0; current_bottom_col_index < bottom->num_columns(); ++current_bottom_col_index) {
            result->set(top->num_rows()+current_bottom_row_index, current_bottom_col_index, bottom->get(current_bottom_row_index, current_bottom_col_index));
        }
    }

    return result;
}

std::shared_ptr<sparse_matrix> sparse_matrix::vertical(const std::shared_ptr<sparse_matrix> bottom)  const {
    assert(num_columns() == bottom->num_columns());
    std::shared_ptr<sparse_matrix> result(new sparse_matrix(num_rows() + bottom->num_rows(), 
                            num_columns(), 
                            num_nnz() + bottom->num_nnz()));
    
    size_t current_nnz_index=0;
    std::vector<double> bottom_vals        = bottom->get_vals();
    std::vector<   int> bottom_row_indices = bottom->get_row_indices();
    std::vector<   int> bottom_col_starts  = bottom->get_col_starts();
    
    // std::vector<double> *result_vals = new std::vector<double>(num_nnz() + bottom->num_nnz());
    // std::vector<   int> *result_row_indices = new std::vector<int>(num_nnz() + bottom->num_nnz());
    // std::vector<   int> *result_col_starts = new std::vector<int>(num_columns()+1);
        
    for (size_t col_start_index=0; col_start_index < num_columns(); ++col_start_index) {
        // std::cout << "col_start_index: " <<col_start_index << std::endl;
        // std::cout << "              a: " << lower_var.col_starts_[col_start_index] << " -> " << lower_var.col_starts_[col_start_index+1] << std::endl;
        // std::cout << "              b: " << upper_var.col_starts_[col_start_index] << " -> " << upper_var.col_starts_[col_start_index+1] << std::endl;
        for (std::size_t current_a_index=col_starts_[col_start_index]; current_a_index < col_starts_[col_start_index+1]; ++current_a_index ) {
            // std::cout << "reading an A" << std::endl;
            result->vals_[current_nnz_index] = vals_[current_a_index];
            result->row_indices_[current_nnz_index] = row_indices_[current_a_index];
            ++current_nnz_index;
        }
        for (std::size_t current_b_index=bottom_col_starts[col_start_index]; current_b_index < bottom_col_starts[col_start_index+1]; ++current_b_index ) {
            // std::cout << "reading a B..." << std::endl;
            result->vals_[current_nnz_index] = bottom_vals[current_b_index];
            result->row_indices_[current_nnz_index] = num_rows() + bottom_row_indices[current_b_index];
            ++current_nnz_index;
        }
        result->col_starts_[col_start_index+1] = current_nnz_index;
    }
    
    // free(vals_);
    // &vals_ = result_vals;
    // free(row_indices_);
    // &row_indices_ = result_row_indices;
    // free(col_starts_);
    // &col_starts_ = result_col_starts;
    
    // std::cout << "result:"<< std::endl << result << std::endl;
    return result;
}

// std::shared_ptr<sparse_matrix> vertical(const std::shared_ptr<sparse_matrix> top, const std::shared_ptr<sparse_matrix> bottom) {
//     assert(top->num_columns() == bottom->num_columns());
//     std::shared_ptr<sparse_matrix> result(new sparse_matrix(top->num_rows() + bottom->num_rows(), 
//                             top->num_columns(), 
//                             top->num_nnz() + bottom->num_nnz()));
//     size_t current_a_index=0, current_b_index=0;
// 
//     size_t current_nnz_index=0;
//     std::vector<double> top_vals        = top->get_vals();
//     std::vector<   int> top_row_indices = top->get_row_indices();
//     std::vector<   int> top_col_starts  = top->get_col_starts();
//     std::vector<double> bottom_vals        = bottom->get_vals();
//     std::vector<   int> bottom_row_indices = bottom->get_row_indices();
//     std::vector<   int> bottom_col_starts  = bottom->get_col_starts();
//     
//     std::vector<double> result_vals(top->num_nnz() + bottom->num_nnz());
//     std::vector<   int> result_row_indices(top->num_nnz() + bottom->num_nnz());
//     std::vector<   int> result_col_starts(top->num_columns()+1);
//         
//     for (size_t col_start_index=0; col_start_index < top->num_columns(); ++col_start_index) {
//         // std::cout << "col_start_index: " <<col_start_index << std::endl;
//         // std::cout << "              a: " << lower_var.col_starts_[col_start_index] << " -> " << lower_var.col_starts_[col_start_index+1] << std::endl;
//         // std::cout << "              b: " << upper_var.col_starts_[col_start_index] << " -> " << upper_var.col_starts_[col_start_index+1] << std::endl;
//         for (std::size_t current_a_index=top_col_starts[col_start_index]; current_a_index < top_col_starts[col_start_index+1]; ++current_a_index ) {
//             // std::cout << "reading an A" << std::endl;
//             result_vals[current_nnz_index] = top_vals[current_a_index];
//             result_row_indices[current_nnz_index] = top_row_indices[current_a_index];
//             ++current_nnz_index;
//         }
//         for (std::size_t current_b_index=bottom_col_starts[col_start_index]; current_b_index < bottom_col_starts[col_start_index+1]; ++current_b_index ) {
//             // std::cout << "reading a B..." << std::endl;
//             result_vals[current_nnz_index] = bottom_vals[current_b_index];
//             result_row_indices[current_nnz_index] = top->num_rows() + result_row_indices[current_b_index];
//             ++current_nnz_index;
//         }
//         result_col_starts[col_start_index+1] = current_nnz_index;
//     }
// 
//     // std::cout << "result:"<< std::endl << result << std::endl;
//     return result;
// }

// TODO move this into dense_matrix
std::shared_ptr<dense_matrix> horizontal(const std::shared_ptr<dense_matrix> left, const std::shared_ptr<dense_matrix> right) {
    assert(left->num_rows() == right->num_rows());
    std::shared_ptr<dense_matrix> result(new dense_matrix(   left->num_rows(),
                            left->num_columns() + right->num_columns()));

    for (size_t current_row_index=0; current_row_index < left->num_rows(); ++current_row_index) {
        for (size_t current_left_col_index=0; current_left_col_index < left->num_columns(); ++current_left_col_index) {
            result->set(current_row_index, current_left_col_index, left->get(current_row_index, current_left_col_index));
        }
        for (size_t current_right_col_index=0; current_right_col_index < right->num_columns(); ++current_right_col_index) {
            result->set(current_row_index, current_right_col_index, right->get(current_row_index, current_right_col_index));
        }
    }

    return result;
}


std::shared_ptr<sparse_matrix> sparse_matrix::horizontal(const std::shared_ptr<sparse_matrix> right) const {
    assert(num_rows() == right->num_rows());
    std::shared_ptr<sparse_matrix> result(new sparse_matrix(   num_rows(),
                            num_columns() + right->num_columns(), 
                            num_nnz() + right->num_nnz()));
    size_t result_nonzero_index=0;
    for (size_t left_nonzeros=0; left_nonzeros < num_nnz(); ++left_nonzeros) {
        result->vals_[result_nonzero_index] = vals_[left_nonzeros];
        result->row_indices_[result_nonzero_index] = row_indices_[left_nonzeros];
        ++result_nonzero_index;
    }
    for (size_t right_nonzeros=0; right_nonzeros < right->num_nnz(); ++right_nonzeros) {
        result->vals_[result_nonzero_index] = right->vals_[right_nonzeros];
        result->row_indices_[result_nonzero_index] = right->row_indices_[right_nonzeros];
        ++result_nonzero_index;
    }
    
    size_t result_column_index=0;
    for (size_t left_cols=0; left_cols < num_columns()+1; ++left_cols) {
        result->col_starts_[result_column_index] = col_starts_[left_cols];
        ++result_column_index;
    }
    result_column_index=num_columns();
    for (size_t right_cols=0; right_cols < right->num_columns()+1; ++right_cols) {
        result->col_starts_[result_column_index] = num_nnz() + right->col_starts_[right_cols];
        ++result_column_index;
    }
    return result;
}

void sparse_matrix::regularize(double hessian_shift, double last_shift) {
    // only allow regularization on square matrices.
    assert(num_rows() == num_columns());
    
    // assert(false); // this is fucked because AMPL doesn't always include the diags.
    
    // std::cout << "trying to regularize: " << hessian_shift << std::endl;
    // for (size_t column_index=0; column_index < num_columns(); ++column_index) {
//             // std::cout << " - col: " << column_index << std::endl;
//             for (size_t nonzero_index=col_starts_[column_index]; nonzero_index < col_starts_[column_index+1]; ++nonzero_index) {
//                 // std::cout << " - col: " << column_index << "; nonzero: " << nonzero_index << "; row:" << row_indices_[nonzero_index] << "; val:" << vals_[nonzero_index] << std::endl;
//                 if (row_indices_[nonzero_index] == column_index) {
//                     vals_[nonzero_index] = vals_[nonzero_index] + hessian_shift - last_shift;
//                     // std::cout << " - col: " << column_index << "; nonzero: " << nonzero_index << "; row:" << row_indices_[nonzero_index] << "; val:" << vals_[nonzero_index] << "UPDATED!" << std::endl;
//                 }
//             }
//         }
//         // last_hessian_shift_ = hessian_shift;
    std::shared_ptr<sparse_matrix> scaled_identity(new sparse_matrix(num_columns(), hessian_shift - last_shift));
    // sparse_matrix sib(num_columns(), hessian_shift - last_shift);
     // *scaled_identity_nosmart = scaled_identity.get();
    
    // scaled_identity->print(std::cout);
    // std::cout << "scaled identity: " << scaled_identity->print(std::cout) << std::endl;
    // std::cout << "scaled identity: " << std::endl << *scaled_identity << std::endl;
    sum(scaled_identity);
    // std::cout << "*this: " << std::endl << *this << std::endl;
    
    // assert(false);
}