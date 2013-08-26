#ifndef SUBPROBLEM_HH_FOZWW1AF
#define SUBPROBLEM_HH_FOZWW1AF

#include <iostream>
#include <vector>

#include "nlp.hh"
#include "iterate.hh"
#include "nlp_state.hh"

class iSQOQuadraticSubproblem : public FunctionWithNLPState {
public:
	iSQOQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate);
	
	void inc_regularization(double hessian_shift);
	
	void print() const;
	int num_qp_variables_, num_qp_constraints_;
	int num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_;
	matrix hessian_;
	bool first_shift_;
	std::vector<double> unshifted_hessian_diagonal_;
	matrix jacobian_;
	std::vector<double> gradient_;
	std::vector<double> lower_bound_;
	std::vector<double> upper_bound_;
	std::vector<double> jacobian_lower_bound_;
	std::vector<double> jacobian_upper_bound_;
	matrix nlp_hessian_;
	matrix nlp_eq_jacobian_;
	matrix nlp_ieq_jacobian_;
	std::vector<double> nlp_objective_gradient_;
private:
protected:
};

#endif /* end of include guard: SUBPROBLEM_HH_FOZWW1AF */
