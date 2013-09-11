
#include "linear_model.hh"
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
        std::cout << "iterate: " << iterate << std::endl;
        std::cout << "nlp eq: " << nlp_->constraints_equality(iterate) << std::endl;
        std::cout << "nlp ieq: " << nlp_->constraints_inequality(iterate) << std::endl;
        
        LinearModelFunction lmf(*nlp_);
        iSQOStep zero_step(nlp_->num_primal(), nlp_->num_dual_eq(), nlp_->num_dual_ieq(), -13, -13);
    	std::cout << " - : " << lmf(iterate, zero_step) << std::endl;
    	std::cout << " -      : " << lmf(iterate, step) << std::endl;
        std::cout << " -      : " << lmf(iterate, zero_step)-lmf(iterate, step) << std::endl;
    	
        
    	std::cout << " - penalty_scaled_objective_decrease: " << penalty_scaled_objective_decrease << std::endl;
    	std::cout << " - iterate_constraint_violation     : " << iterate_constraint_violation << std::endl;
    	std::cout << " - iterate_step_constraint_violation: " << iterate_step_constraint_violation << std::endl;
    	std::cout << " - retval                           : " << retval << std::endl;
        
        printf("psod: %9.2e\ticv: %9.2e\tiscv: %9.2e\t\t\tlinear reduction%9.2e\n", penalty_scaled_objective_decrease, iterate_constraint_violation, iterate_step_constraint_violation, retval);
    }
	return retval;
}
