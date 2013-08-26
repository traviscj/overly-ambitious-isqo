
#include "linear_model_reduction.hh"

LinearReductionFunction::LinearReductionFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
	
}
double LinearReductionFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
	bool PRINT=false;
	// double f = nlp_->objective(iterate);
	std::vector<double> gradient = nlp_->objective_gradient(iterate);
	// dotproduct(gradient, )
	double dot_product=0.0;
	for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index) {
		dot_product += gradient[primal_index]*step.primal_values_[primal_index];
	}
	
	if (PRINT) std::cout << "linear decrease: " << std::endl;
	if (PRINT) std::cout << " - " << -iterate.penalty_parameter_*dot_product << std::endl;
	if (PRINT) std::cout << " - " << constraint_violation_func_(iterate) << std::endl;
	if (PRINT) std::cout << " - " << constraint_violation_func_(iterate,step) << std::endl;
	if (PRINT) std::cout << " - " << -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step) << std::endl;
	return -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step);
}