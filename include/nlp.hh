
#ifndef __GUARD_NLP_HH
#define __GUARD_NLP_HH

#include <memory>

#include "iterate.hh"
#include "matrix.hh"
#include "utilities.hh"

class Nlp {
	// this class implements 
public:
	Nlp(size_t num_primal, size_t num_dual_eq, size_t num_dual_ieq);
	
	// starting point:
	virtual iSQOIterate initial() = 0;
	
	// zeroth order NLP quantities:
	virtual double objective(const iSQOIterate &iterate) = 0;
	virtual std::vector<double> constraints_equality(const iSQOIterate &iterate) = 0;
	virtual std::vector<double> constraints_inequality(const iSQOIterate &iterate) = 0;
	
	// first order NLP quantities:
	virtual std::vector<double> objective_gradient(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate) = 0;
	virtual std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate) = 0;
    // virtual sparse_matrix constraints_equality_jacobian_sparse(const iSQOIterate &iterate) = 0;
    // virtual sparse_matrix constraints_inequality_jacobian_sparse(const iSQOIterate &iterate) = 0;
    
	// second order NLP quantities:
	virtual std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate) = 0;
    // virtual sparse_matrix lagrangian_hessian_sparse(const iSQOIterate &iterate) = 0;

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
