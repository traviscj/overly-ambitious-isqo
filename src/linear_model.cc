
#include "linear_model.hh"
#include "utilities.hh"

LinearModelFunction::LinearModelFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
	
}
double LinearModelFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) {
	double f = nlp_->objective(iterate);
	std::vector<double> gradient = nlp_->objective_gradient(iterate);
	// dotproduct(gradient, )
	double dot_product=0.0;
	for (size_t primal_index=0; primal_index < nlp_->num_primal(); ++primal_index) {
		dot_product += gradient[primal_index]*step.get_primal_value(primal_index);
	}
	return iterate.get_penalty_parameter()*(f + dot_product) + constraint_violation_func_(iterate, step);
}