
#ifndef __GUARD_MATRIX_HH
#define __GUARD_MATRIX_HH

#include <iostream>
#include <memory>
#include <vector>
#include <cassert>

#include "utilities.hh"

//! \brief a matrix
//!
//! 
class matrix_base_class {
public:
    matrix_base_class(std::size_t rows, std::size_t columns);
    virtual ~matrix_base_class() = 0;

	virtual std::ostream &print(std::ostream& os) const = 0;
	virtual std::vector<double> multiply(const std::vector<double> &x) const = 0;
	virtual std::vector<double> multiply_transpose(const std::vector<double> &y) const = 0;
    virtual void regularize(double hessian_shift, double last_shift) = 0;
    
    std::size_t num_columns() const { return columns_;}
    std::size_t num_rows() const { return rows_;}
    
private:
protected:
    std::size_t rows_;
    std::size_t columns_;
    double last_hessian_shift_;
};

//! \brief a matrix
//!
//! 
class dense_matrix : public matrix_base_class {
public:
    // dense_matrix(std::size_t rows, std::size_t columns);
    
    //! \brief create an empty dense matrix
    //!
    //! \param num_rows the number of rows the 'null' sparse matrix has
    //! \param num_cols the number of columns the 'null' sparse matrix has
    dense_matrix(std::size_t num_rows, std::size_t num_cols);
    
    // 
    //! \brief a 'relevant identity' dense matrix constructor
    //!
    //! \param num_cols the number of columns the sparse matrix has
    //! \param relevant_variables specifies which rows of the identity to include.
    //! \param cur_sign value to assign the nonzero elements
    dense_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign);
    
    //! \brief pure identity dense matrix constructor
    //! \param num_variables size of the matrix to be created
    //! \param scalar value to assign along diagonal
    //!
    //! This will create an identity matrix, of size num_variables, with the sparse structure:
    //!  - values_{I_n}: scalar*[1.0, 1.0, 1.0,...] (num_variables entries)
    //!  - row_indices_{I_n}: [0,1,2,...,num_variables-1] (num_variables entries)
    //!  - col_starts_{I_n}: [0,1,2,...,num_variables] (num_variables+1 entries)
    dense_matrix(std::size_t num_variables, double scalar);
    
    void regularize(double hessian_shift, double last_shift) {
        // only allow regularization on square matrices.
        assert(num_rows() == num_columns());
        for (size_t diagonal_entry=0; diagonal_entry<num_rows(); ++diagonal_entry) {
                set(diagonal_entry,diagonal_entry, (hessian_shift - last_shift) + get(diagonal_entry, diagonal_entry));
        }
    }
    void set(std::size_t r, std::size_t c, double val);
	double get(std::size_t r, std::size_t c) const;

	std::ostream &print(std::ostream& os) const;
	std::vector<double> multiply(const std::vector<double> &x) const;
	std::vector<double> multiply_transpose(const std::vector<double> &y) const;

	std::vector<double> data_;
private:
protected:
	
};

inline std::ostream& operator<< (std::ostream& os, const dense_matrix& m) {
	m.print(os);
    return os;
}


class sparse_matrix : public matrix_base_class {
public:
    //! \brief create a 'null' sparse matrix
    //!
    //! \param num_rows the number of rows the 'null' sparse matrix has
    //! \param num_cols the number of columns the 'null' sparse matrix has
    //! \param num_nonzeros the number of nonzeros the matrix will have.
    sparse_matrix(std::size_t num_rows, std::size_t num_cols, std::size_t num_nonzeros);
    
    // 
    //! \brief a 'relevant identity' sparse matrix constructor
    //!
    //! \param num_cols the number of columns the sparse matrix has
    //! \param relevant_variables specifies which rows of the identity to include.
    //! \param cur_sign value to assign the nonzero elements
    sparse_matrix(std::size_t num_cols, std::vector<std::size_t> relevant_variables, double cur_sign);
    
    //! \brief pure identity sparse matrix constructor
    //! \param num_variables size of the matrix to be created
    //! \param scalar value to assign along diagonal
    //!
    //! This will create an identity matrix, of size num_variables, with the sparse structure:
    //!  - values_{I_n}: scalar*[1.0, 1.0, 1.0,...] (num_variables entries)
    //!  - row_indices_{I_n}: [0,1,2,...,num_variables-1] (num_variables entries)
    //!  - col_starts_{I_n}: [0,1,2,...,num_variables] (num_variables+1 entries)
    sparse_matrix(std::size_t num_variables, double scalar);
    
    std::vector<double> multiply(const std::vector<double> &) const;
    std::vector<double> multiply_transpose(const std::vector<double> &) const;
    std::ostream &print(std::ostream& os) const;
    
    void regularize(double hessian_shift, double last_shift) {
        // only allow regularization on square matrices.
        assert(num_rows() == num_columns());
        
        // std::cout << "trying to regularize: " << hessian_shift << std::endl;
        for (size_t column_index=0; column_index < num_columns(); ++column_index) {
            // std::cout << " - col: " << column_index << std::endl;
            for (size_t nonzero_index=col_starts_[column_index]; nonzero_index < col_starts_[column_index+1]; ++nonzero_index) {
                // std::cout << " - col: " << column_index << "; nonzero: " << nonzero_index << "; row:" << row_indices_[nonzero_index] << "; val:" << vals_[nonzero_index] << std::endl;
                if (row_indices_[nonzero_index] == column_index) {
                    vals_[nonzero_index] = vals_[nonzero_index] + hessian_shift - last_shift;
                    // std::cout << " - col: " << column_index << "; nonzero: " << nonzero_index << "; row:" << row_indices_[nonzero_index] << "; val:" << vals_[nonzero_index] << "UPDATED!" << std::endl;
                }
            }
        }
        // last_hessian_shift_ = hessian_shift;
    }
    
    std::size_t num_nnz() const { return vals_.size();}

    std::vector<double> vals_;
    std::vector<int> row_indices_;
    std::vector<int> col_starts_;
};

//! \brief print a sparse matrix
//!
//! here, we print out the underlying vectors associated with the sparse matrix 'mat', to the output stream 'os'.
inline std::ostream& operator<< (std::ostream& os, const sparse_matrix& mat) {
    return mat.print(os);
}

//! \brief print a sparse matrix
//!
//! here, we print out the underlying vectors associated with the sparse matrix 'mat', to the output stream 'os'.
inline std::ostream& operator<< (std::ostream& os, const std::shared_ptr<matrix_base_class>& mat) {
    return mat->print(os);
}


// TODO: 'add_to_diag' function for sparse_matrix.
// TODO: change this to 'append_below' & put into sparse_matrix class.
std::shared_ptr<dense_matrix> vertical(const std::shared_ptr<dense_matrix> top, const std::shared_ptr<dense_matrix> bottom);
std::shared_ptr<sparse_matrix> vertical(const std::shared_ptr<sparse_matrix> top, const std::shared_ptr<sparse_matrix> bottom);

// TODO: change this to 'append_right' & put into sparse_matrix class.
std::shared_ptr<dense_matrix> horizontal(const std::shared_ptr<dense_matrix> left, const std::shared_ptr<dense_matrix> right);
std::shared_ptr<sparse_matrix> horizontal(const std::shared_ptr<sparse_matrix> left, const std::shared_ptr<sparse_matrix> right);

#endif
