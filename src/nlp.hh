
#ifndef __GUARD_NLP_HH
#define __GUARD_NLP_HH

#include <memory>

#include "iterate.hh"
#include "matrix.hh"
#include "utilities.hh"

//! \brief class returning NLP evaluations.
//!
//! We assume that the Nlp object returns a model of the form (for a specified \f$x\f$, \f$\lambda_{\mathcal{E}}\f$, and \f$\lambda_{\mathcal{I}}\f$)
//! \f{align*}{\min& \qquad  f(x) + \nabla f(x)^T d + \frac{1}{2}d^T\left[ \nabla^2_{xx} f(x) + \sum_{i\in\mathcal{E}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} + \sum_{i\in\mathcal{I}} \nabla^2_{xx} c_{i}(x) \cdot \lambda_{(i)} \right]d \\ \text{s.t.}& \qquad  J_{\mathcal{E}}\cdot d + c_{\mathcal{E}}(x) = 0\\& \qquad   J_{\mathcal{I}}\cdot d + c_{\mathcal{I}}(x) \leq 0 \f}

class Nlp {
	// this class implements 
public:
	Nlp(size_t num_primal, size_t num_dual_eq, size_t num_dual_ieq);
	
	// starting point:
	virtual iSQOIterate initial(double penalty_parameter) = 0;
	
	// zeroth order NLP quantities:
	virtual double objective(const iSQOIterate &iterate) = 0;
	virtual std::vector<double> constraints_equality(const iSQOIterate &iterate) = 0;
	virtual std::vector<double> constraints_inequality(const iSQOIterate &iterate) = 0;
	
	// first order NLP quantities:
	virtual std::vector<double> objective_gradient(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate) = 0;
    
	// second order NLP quantities:
	virtual std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate) = 0;

	size_t num_primal();
	size_t num_dual();
	size_t num_dual_eq();
	size_t num_dual_ieq();
protected:
	size_t num_primal_;
	size_t num_dual_eq_;
	size_t num_dual_ieq_;
};

#endif
