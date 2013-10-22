
#include <iostream>
#include <memory>
using namespace std;

#include "matrix.hh"



void build_matrix(std::size_t n_eq, std::size_t n_ieq) {
    std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol0(sparse_matrix::diagonal_matrix( n_eq, n_eq, -1.0));
    std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol1(new sparse_matrix( n_eq, n_ieq, 0));
    std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol2(sparse_matrix::diagonal_matrix( n_eq, n_eq,  +1.0) );
    
    std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol0(new sparse_matrix(n_ieq,  n_eq,      0));
    std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol1(sparse_matrix::diagonal_matrix(n_ieq, n_ieq, -1.0));
    std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol2(new sparse_matrix(n_ieq,  n_eq,      0));
    
    // cout << jacobian_slack_brow0_bcol0 << endl;
    // cout << jacobian_slack_brow0_bcol1 << endl;
    // cout << jacobian_slack_brow0_bcol2 << endl;
    // cout << jacobian_slack_brow1_bcol0 << endl;
    // cout << jacobian_slack_brow1_bcol1 << endl;
    // cout << jacobian_slack_brow1_bcol2 << endl;
    
    // // // build the all-rows slack matrices:
    std::shared_ptr<sparse_matrix> jacobian_slacks_column0 = jacobian_slack_brow0_bcol0->vertical(jacobian_slack_brow1_bcol0);
    std::shared_ptr<sparse_matrix> jacobian_slacks_column1 = jacobian_slack_brow0_bcol1->vertical(jacobian_slack_brow1_bcol1);
    std::shared_ptr<sparse_matrix> jacobian_slacks_column2 = jacobian_slack_brow0_bcol2->vertical(jacobian_slack_brow1_bcol2);
    // 
    // cout << "jacobian_slacks_column0" << endl << jacobian_slacks_column0 << endl;
    // cout << "jacobian_slacks_column1" << endl << jacobian_slacks_column1 << endl;
    // cout << "jacobian_slacks_column2" << endl << jacobian_slacks_column2 << endl;
    
    
    // 
    // horizontally-concatenate columns to make the full slack matrix:
    std::shared_ptr<sparse_matrix> jacobian_slacks_column0_column1 = jacobian_slacks_column0->horizontal(jacobian_slacks_column1);
    std::shared_ptr<sparse_matrix> jacobian_slacks_all             = jacobian_slacks_column0_column1->horizontal(jacobian_slacks_column2);

    // cout << "jacobian_slacks_column0_column1" << endl << jacobian_slacks_column0_column1 << endl;
    // cout << "jacobian_slacks_all" << endl 
    cout << "mat = zeros(" << n_eq + n_ieq << ", " << 2*n_eq + n_ieq << ");" << endl;
    cout << jacobian_slacks_all << endl;
}
int main() {
    cout << "Here's the matrix test!" << endl;

    std::size_t n_eq = 3, n_ieq=5;
    
    // static std::shared_ptr<dense_matrix> diagonal_matrix(std::size_t rows, std::size_t columns, double diagonal, std::vector<int> *selected_rows=NULL);
    
    
    build_matrix(0,0);
    
    build_matrix(0,n_ieq);
    
    build_matrix(n_eq,0);
    
    build_matrix(n_eq,n_ieq);
    
    cout << " ------------- next phase" << endl;
    std::shared_ptr<sparse_matrix> no_row_test(new sparse_matrix( 0, n_ieq, 0));
    std::shared_ptr<sparse_matrix> no_col_test(new sparse_matrix( n_eq, 0, 0));
    
    cout << no_row_test << endl;
    cout << no_col_test << endl;
    
    return 0;
}