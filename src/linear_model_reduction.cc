
#include "linear_model_reduction.hh"

LinearReductionFunction::LinearReductionFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
	
}
double LinearReductionFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
	bool PRINT=false;
	std::vector<double> gradient = nlp_->objective_gradient(iterate);
    double penalty_scaled_objective_decrease = -iterate.get_penalty_parameter()*dot_product(gradient, step.get_primal_values());
    double iterate_constraint_violation = constraint_violation_func_(iterate);
    double iterate_step_constraint_violation = constraint_violation_func_(iterate,step);
    double retval =  penalty_scaled_objective_decrease + iterate_constraint_violation - iterate_step_constraint_violation;
    if (PRINT && (retval < 0)) {
    	std::cout << std::endl;
        std::cout << "linear decrease was negative!!" << std::endl;
    	std::cout << " - penalty_scaled_objective_decrease: " << penalty_scaled_objective_decrease << std::endl;
    	std::cout << " - iterate_constraint_violation     : " << iterate_constraint_violation << std::endl;
    	std::cout << " - iterate_step_constraint_violation: " << iterate_step_constraint_violation << std::endl;
    	std::cout << " - retval                           : " << retval << std::endl;
    }
	return retval;
}
