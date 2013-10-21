#include <memory>
#include <cmath>

#include "utilities.hh"
#include "residual_function.hh"

ResidualFunction::ResidualFunction(iSQOControlPanel &control, Nlp &nlp) : FunctionWithNLPState(control, nlp) {
	// std::cout << "-- Initializing nlp residual function." << std::endl;
}
double ResidualFunction::resid_helper(const iSQOIterate &iterate, const std::vector<double> &stationarity, const std::vector<double> &constraint_eq_values, const std::vector<double> &constraint_ieq_values, const std::vector<double> &constraint_eq_dual_values, const std::vector<double> &constraint_ieq_dual_values) const {
	bool PRINT=true;
    // std::cout << "******* resid_helper" << std::endl;
	std::vector<double> rho(stationarity.size() + 2*constraint_eq_values.size() + 2*constraint_ieq_values.size());
	
	std::shared_ptr<matrix_base_class> jac_ce = nlp_->constraints_equality_jacobian(iterate);
	std::shared_ptr<matrix_base_class> jac_ci = nlp_->constraints_inequality_jacobian(iterate);
			
	size_t rho_index = 0;

	double eq_signflip = 1;
	double ieq_signflip = 1;
    if (PRINT) std::cout << "                                                                              sta pre: " << stationarity << std::endl;
	std::vector<double> eq_jacobian_trans_times_mults = jac_ce->multiply_transpose(constraint_eq_dual_values);
    if (PRINT) std::cout << "eq_jacobian' .constraint_eq_dual_values" << constraint_eq_dual_values << ":                                            " << eq_jacobian_trans_times_mults << std::endl;
    
	std::vector<double> ieq_jacobian_trans_times_mults = jac_ci->multiply_transpose(constraint_ieq_dual_values);
    if (PRINT) std::cout << "ieq_jacobian'.constraint_ieq_dual_values" << constraint_ieq_dual_values << ": " << ieq_jacobian_trans_times_mults << std::endl;
	
	for (size_t stationarity_index=0; stationarity_index < stationarity.size(); ++stationarity_index) {
		rho[rho_index] = stationarity[stationarity_index] + eq_jacobian_trans_times_mults[stationarity_index] + ieq_jacobian_trans_times_mults[stationarity_index];
        if (PRINT) std::cout << "******** Stationarity, i=" << stationarity_index << ": " << rho[rho_index] << std::endl;
		++rho_index;
	}
    if (PRINT) std::cout << "                                                                              sta pre: " << stationarity << std::endl;
    
    std::cout << " constraint_eq_values: " << constraint_eq_values  << "\t" <<  constraint_eq_dual_values << std::endl;
    std::cout << "constraint_ieq_values: " << constraint_ieq_values << "\t" << constraint_ieq_dual_values << std::endl;
	for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {        
		rho[rho_index] = fmin(bracket_plus(constraint_eq_values[constraint_eq_index]), 1.0 - eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
		++rho_index;
    }
	for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {        
		rho[rho_index] = fmin(bracket_minus(constraint_eq_values[constraint_eq_index]), 1.0 + eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
		++rho_index;
	}
    // if (PRINT) std::cout << "constraint_ieq_dual_values: " << constraint_ieq_dual_values << std::endl;
	for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
		rho[rho_index] = fmin(bracket_plus(constraint_ieq_values[constraint_ieq_index]), 1.0 - ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
		++rho_index;
    }
    for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
		rho[rho_index] = fmin(bracket_minus(constraint_ieq_values[constraint_ieq_index]), 0.0 + ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
		++rho_index;
		
	}
    if (PRINT) std::cout << "******** Full Rho: " << rho << std::endl;
    if (PRINT) std::cout << "rho(" << iterate.get_penalty_parameter() << "): " << rho << std::endl;
	double total = 0.0;
	for (size_t rho_index = 0; rho_index < rho.size(); ++rho_index)
		total += rho[rho_index]*rho[rho_index];
    if (PRINT) std::cout << "******** ||Rho||_2: " << sqrt(total) << std::endl;
	
	return sqrt(total);
}
double ResidualFunction::operator()(const iSQOIterate &iterate) const {
    // std::cout << "******* ResidualFunction::operator()(const iSQOIterate &iterate) const {" << std::endl;
	bool PRINT=false;
    if (PRINT) std::cout << "OPERATOR FOR RESID FUNC: " << iterate.get_penalty_parameter() << std::endl;
	
	std::vector<double> stationarity = nlp_->objective_gradient(iterate);
	
	for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
		// mu*grad f:
		stationarity[stationarity_index] *= iterate.get_penalty_parameter();
	}
	
	std::vector<double> eq_con_values = nlp_->constraints_equality(iterate);
	std::vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
	
	if (PRINT) {
		std::cout << std::endl << "mu: " << iterate.get_penalty_parameter() << std::endl;
		std::cout << "sta (pre): " << nlp_->objective_gradient(iterate) << std::endl;
		std::cout << "sta: " << stationarity << std::endl;
		std::cout << "constraints  eq: " << eq_con_values << std::endl;
		std::cout << "constraints ieq: " << ieq_con_values << std::endl;
	
		std::cout << "jac eq: " << nlp_->constraints_equality_jacobian(iterate) << std::endl;
	
		std::cout << "       dual  eq: " << iterate.get_dual_eq_values() << std::endl;
		std::cout << "       dual ieq: " << iterate.get_dual_ieq_values() << std::endl;
	}
    std::vector<double> dual_eq_values(nlp_->num_dual_eq());
    dual_eq_values.assign(iterate.get_dual_eq_values()->begin(), iterate.get_dual_eq_values()->end());
    std::vector<double> dual_ieq_values(nlp_->num_dual_ieq());
    dual_ieq_values.assign(iterate.get_dual_ieq_values()->begin(), iterate.get_dual_ieq_values()->end());
    
	return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, dual_eq_values, dual_ieq_values);
}
double ResidualFunction::operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
    // std::cout << "******* ResidualFunction::operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step" << std::endl;
	bool PRINT=true;
    
	if (PRINT) std::cout << "===== START OF STEP RESIDUAL CALCULATION" << std::endl;

    // std::cout << "step: " << step << std::endl;

	std::vector<double> stationarity = nlp_->objective_gradient(iterate);	
	std::shared_ptr<matrix_base_class> hessian = nlp_->lagrangian_hessian(iterate);
	std::vector<double> hessian_times_step = hessian->multiply(step.get_primal_values());
	
    for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
		// mu*grad f + H*d:
		stationarity[stationarity_index] = iterate.get_penalty_parameter()*(stationarity[stationarity_index]) + hessian_times_step[stationarity_index];
		
	}
	
	std::shared_ptr<matrix_base_class> jacobian_eq = nlp_->constraints_equality_jacobian(iterate);
	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
    if (PRINT) std::cout << "con_values_eq 0: " << con_values_eq << std::endl;
    if (PRINT) std::cout << "jacobian_eq: " << jacobian_eq << std::endl;
    if (PRINT) std::cout << "step.get_primal_values(): " << step.get_primal_values() << std::endl;
	// std::cout << "about to do j eq" << std::endl;
	std::vector<double> linear_step_eq = jacobian_eq->multiply(step.get_primal_values());
    if (PRINT) std::cout << "linear_step_eq 0: " << linear_step_eq << std::endl;
	
	// std::cout << "about to do j ieq" << std::endl;
	std::shared_ptr<matrix_base_class> jacobian_ieq = nlp_->constraints_inequality_jacobian(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
    if (PRINT) std::cout << "con_values_ieq 0: " << con_values_ieq << std::endl;
    if (PRINT) std::cout << "jacobian_ieq: " << jacobian_ieq << std::endl;
    if (PRINT) std::cout << "step.get_primal_values(): " << step.get_primal_values() << std::endl;
	std::vector<double> linear_step_ieq = jacobian_ieq->multiply(step.get_primal_values());
    if (PRINT) std::cout << "linear_step_ieq 0: " << linear_step_ieq << std::endl;

	for (size_t eq_con_index=0; eq_con_index< nlp_->num_dual_eq(); ++eq_con_index) {
        std::cout << "equality copy. con_values_eq[eq_con_index=" << eq_con_index << "]: " << con_values_eq[eq_con_index] << " + "<< linear_step_eq[eq_con_index];
		con_values_eq[eq_con_index] += linear_step_eq[eq_con_index];
        std::cout << " = " << con_values_eq[eq_con_index] << std::endl;
	}
	for (size_t ieq_con_index=0; ieq_con_index< nlp_->num_dual_ieq(); ++ieq_con_index) {
        std::cout << "inequality copy. con_values_ieq[ieq_con_index=" << ieq_con_index << "]: " << con_values_ieq[ieq_con_index] << " + "<< linear_step_ieq[ieq_con_index];
		con_values_ieq[ieq_con_index] += linear_step_ieq[ieq_con_index];
        std::cout << " = " << con_values_ieq[ieq_con_index] << std::endl;
	}
    if (PRINT) std::cout << "con_values_ieq 1: " << con_values_ieq << std::endl;
	if (PRINT) {
        std::cout << "step: " << step << std::endl;
		std::cout << std::endl << "mu: " << iterate.get_penalty_parameter() << std::endl;
		std::cout << "sta (pre): " << nlp_->objective_gradient(iterate) << std::endl;
		std::cout << "sta: " << stationarity << std::endl;
		std::cout << "constraints  eq: " << con_values_eq << std::endl;
		std::cout << "constraints ieq: " << con_values_ieq << std::endl;
	
        std::cout << "jac  eq: " << nlp_->constraints_equality_jacobian(iterate) << std::endl;
        
        std::cout << "jac ieq: " << nlp_->constraints_inequality_jacobian(iterate) << std::endl;
	
		std::cout << "       dual  eq: " << (iterate.get_dual_eq_values()) << std::endl;
		std::cout << "       dual ieq: " << (iterate.get_dual_ieq_values()) << std::endl;
	}
    std::vector<double> dual_eq_values(nlp_->num_dual_eq());
    dual_eq_values.assign(step.get_dual_eq_values().begin(), step.get_dual_eq_values().end());
    std::vector<double> dual_ieq_values(nlp_->num_dual_ieq());
    dual_ieq_values.assign(step.get_dual_ieq_values().begin(), step.get_dual_ieq_values().end());
    
    double retval = resid_helper(iterate, stationarity, con_values_eq, con_values_ieq, dual_eq_values, dual_ieq_values);
    if (retval > 1e-5) {
        std::cerr << "WHOA THERE, SLOW DOWN, BIG RESIDUAL: " << retval << std::endl;
    }
    if (PRINT) std::cout << "===== END OF STEP RESIDUAL CALCULATION" << std::endl;
	return retval;
}
