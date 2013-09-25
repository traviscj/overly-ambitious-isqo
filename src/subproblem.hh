#ifndef SUBPROBLEM_HH_FOZWW1AF
#define SUBPROBLEM_HH_FOZWW1AF

#include <iostream>
#include <map>
#include <vector>

#include "nlp.hh"
#include "iterate.hh"
#include "nlp_state.hh"
#include "utilities.hh"

// TODO reformulate ieq constraints into bounds on slacks instead of this hacky-ass setup?
// 
//! \brief a class for actually constructing the qpOASES subproblem.
//! 
//! Here, we convert Nlp's format(for a particular \f$x\f$, \f$\lambda_{\mathcal{E}}\f$, and \f$\lambda_{\mathcal{I}}\f$)
//!  \f{align*}{\min_{d}& \qquad  f(x) + \nabla f(x)^T d + \frac{1}{2}d^T\left[ \nabla^2_{xx} f(x) + \sum_{i\in\mathcal{E}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} + \sum_{i\in\mathcal{I}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} \right]d \\ \text{s.t.}& \qquad  J_{\mathcal{E}}\cdot d + c_{\mathcal{E}}(x) = 0\\& \qquad   J_{\mathcal{I}}\cdot d + c_{\mathcal{I}}(x) \leq 0 \f}
//! into
//! \f{align*}{\min_{d}& \qquad \mu f(x) + \mu\nabla f(x)^T d + e^T(p_{\mathcal{E}}+n_{\mathcal{E}}) + e^T(p_{\mathcal{I}}) + \frac{1}{2}d^T\left[ \mu\nabla^2_{xx} f(x) + \sum_{i\in\mathcal{E}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} + \sum_{i\in\mathcal{I}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} \right]d \\ \text{s.t.}& \qquad  lb <= J_{\texttt{qpOASES}}\cdot d - Ip + In <= ub\\& \qquad (p,n) \geq 0 \f}
//! by stacking
//! \f[ J_{\texttt{qpOASES}} = \begin{pmatrix} J_{\mathcal{E}} \\ J_{\mathcal{I}} \end{pmatrix}, lb = \begin{pmatrix} -c_{\mathcal{E}}(x) \\ -\infty \end{pmatrix}, ub = \begin{pmatrix} -c_{\mathcal{E}}(x) \\ -c_{\mathcal{I}}(x) \end{pmatrix} \f]
//! \todo change this to bounds on the slacks instead.

class iSQOQuadraticSubproblem : public FunctionWithNLPState {
public:
	iSQOQuadraticSubproblem(iSQOControlPanel &control, Nlp &nlp, const iSQOIterate &iterate);
    ~iSQOQuadraticSubproblem() {
        // std::cout << "I had " << prior_constructed_hessians_.size() << " items in the hessian cache.\n"
            // << " & " << prior_constructed_jacobians_.size() << " items in the jacobian cache.\n" << std::endl;
    }
	void inc_regularization(double hessian_shift, double last_shift);
    
    size_t num_primal() { return num_qp_variables_; }
    size_t num_dual() { return num_qp_constraints_; }
    
	std::ostream &print(std::ostream &os) const;
	size_t num_qp_variables_, num_qp_constraints_;
	size_t num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_;
	std::shared_ptr<matrix_base_class> hessian_;
    // bool first_shift_;
	std::vector<double> unshifted_hessian_diagonal_;
	std::shared_ptr<matrix_base_class> jacobian_;
	std::vector<double> gradient_;
	std::vector<double> lower_bound_;
	std::vector<double> upper_bound_;
	std::vector<double> jacobian_lower_bound_;
	std::vector<double> jacobian_upper_bound_;
	std::shared_ptr<matrix_base_class> nlp_hessian_;
	std::shared_ptr<matrix_base_class> nlp_eq_jacobian_;
	std::shared_ptr<matrix_base_class> nlp_ieq_jacobian_;
	std::vector<double> nlp_objective_gradient_;

    const iSQOIterate *iterate_pointer;
	
private:
protected:
    void setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<matrix_base_class> nlp_eq_jacobian, std::shared_ptr<matrix_base_class> nlp_ieq_jacobian, std::shared_ptr<matrix_base_class> nlp_hessian);
	
    void setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<dense_matrix> nlp_eq_jacobian, std::shared_ptr<dense_matrix> nlp_ieq_jacobian, std::shared_ptr<dense_matrix> nlp_hessian);
    void setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<sparse_matrix> nlp_eq_jacobian, std::shared_ptr<sparse_matrix> nlp_ieq_jacobian, std::shared_ptr<sparse_matrix> nlp_hessian);
    

    std::map<int, std::shared_ptr<matrix_base_class> > prior_constructed_hessians_;
    std::map<int, std::shared_ptr<matrix_base_class> > prior_constructed_jacobians_;
    // static std::map<int, std::shared_ptr<matrix_base_class> > prior_constructed_hessians;
};

//! \brief print a subproblem
//!
//! here, we print out the underlying vectors associated with the subproblem subproblem, to the output stream 'os'.
inline std::ostream& operator<< (std::ostream& os, const iSQOQuadraticSubproblem& subproblem) {
    return subproblem.print(os);
}

// class iSQOSparseQuadraticSubproblem : public iSQOQuadraticSubproblem {
// public:
//     iSQOSparseQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate);
//     void inc_regularization(double hessian_shift);
//     std::shared_ptr<matrix_base_class> jacobian_sparse_;
//     std::shared_ptr<matrix_base_class> hessian_sparse_;
//     std::shared_ptr<matrix_base_class> nlp_hessian_sparse_;
//     std::shared_ptr<matrix_base_class> nlp_eq_jacobian_sparse_;
//     std::shared_ptr<matrix_base_class> nlp_ieq_jacobian_sparse_;
//     
//     double hessian_shift_;
// private:
// protected:
//     void setup_matrix_data_sparse(const iSQOIterate &);
// };

#endif /* end of include guard: SUBPROBLEM_HH_FOZWW1AF */
