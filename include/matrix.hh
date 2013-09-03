
#ifndef __GUARD_MATRIX_HH
#define __GUARD_MATRIX_HH

#include <iostream>
#include <vector>

#include "utilities.hh"

//! \brief a matrix
//!
//! 
class matrix {
public:
	matrix(std::size_t rows, std::size_t columns);
	void set(std::size_t r, std::size_t c, double val);
	double get(std::size_t r, std::size_t c) const;

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

class sparse_matrix {
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
    
    std::vector<double> multiply(const std::vector<double> &);
    std::vector<double> multiply_transpose(const std::vector<double> &);

    
    std::size_t num_columns() const { return columns_;}
    std::size_t num_rows() const { return rows_;}
    std::size_t num_nnz() const { return vals_.size();}
    std::size_t rows_;
    std::size_t columns_;
    std::vector<double> vals_;
    std::vector<int> row_indices_;
    std::vector<int> col_starts_;
};

//! \brief print a sparse matrix
//!
//! here, we print out the underlying vectors associated with the sparse matrix 'mat', to the output stream 'os'.
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

#endif
