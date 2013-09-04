#include <memory>
#include <cmath>

#include "utilities.hh"
#include "residual_function.hh"

ResidualFunction::ResidualFunction(Nlp &nlp) : FunctionWithNLPState(nlp) {
	// std::cout << "-- Initializing nlp residual function." << std::endl;
}
double ResidualFunction::resid_helper(const iSQOIterate &iterate, std::vector<double> stationarity, std::vector<double> constraint_eq_values, std::vector<double> constraint_ieq_values, std::vector<double> constraint_eq_dual_values, std::vector<double> constraint_ieq_dual_values) const {
	std::vector<double> rho(stationarity.size() + 2*constraint_eq_values.size() + 2*constraint_ieq_values.size());
	
	std::shared_ptr<matrix_base_class> jac_ce = nlp_->constraints_equality_jacobian(iterate);
	std::shared_ptr<matrix_base_class> jac_ci = nlp_->constraints_inequality_jacobian(iterate);
			
	size_t rho_index = 0;

	double eq_signflip = 1;
	double ieq_signflip = 1;
	// std::cout << "sta pre: " << stationarity << std::endl << std::endl << std::endl;
	std::vector<double> eq_jacobian_trans_times_mults = jac_ce->multiply_transpose(constraint_eq_dual_values);
	std::vector<double> ieq_jacobian_trans_times_mults = jac_ci->multiply_transpose(constraint_ieq_dual_values);
	
	for (size_t stationarity_index=0; stationarity_index < stationarity.size(); ++stationarity_index) {
		rho[rho_index] = stationarity[stationarity_index] + eq_jacobian_trans_times_mults[stationarity_index] + ieq_jacobian_trans_times_mults[stationarity_index];
		++rho_index;
	}
	for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {
		rho[rho_index] = std::min(bracket_plus(constraint_eq_values[constraint_eq_index]), 1.0 - eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
		++rho_index;
		
		rho[rho_index] = std::min(bracket_minus(constraint_eq_values[constraint_eq_index]), 1.0 + eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
		++rho_index;
	}
	
	for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
		rho[rho_index] = std::min(bracket_plus(constraint_ieq_values[constraint_ieq_index]), 1.0 - ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
		++rho_index;

		rho[rho_index] = std::min(bracket_minus(constraint_ieq_values[constraint_ieq_index]), 0.0 + ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
		++rho_index;
		
	}
	// std::cout << "rho(" << iterate.penalty_parameter_ << "): " << rho << std::endl;
	double total = 0.0;
	for (size_t rho_index = 0; rho_index < rho.size(); ++rho_index)
		total += rho[rho_index]*rho[rho_index];
	
	return sqrt(total);
}
double ResidualFunction::operator()(const iSQOIterate &iterate) const {
	bool PRINT=false;
	// std::cout << "OPERATOR FOR RESID FUNC: " << iterate.penalty_parameter_ << std::endl;
	
	std::vector<double> stationarity = nlp_->objective_gradient(iterate);
	
	for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
		// mu*grad f:
		stationarity[stationarity_index] *= iterate.penalty_parameter_;
	}
	
	std::vector<double> eq_con_values = nlp_->constraints_equality(iterate);
	std::vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
	
	if (PRINT) {
		std::cout << std::endl << "mu: " << iterate.penalty_parameter_ << std::endl;
		std::cout << "sta (pre): " << nlp_->objective_gradient(iterate) << std::endl;
		std::cout << "sta: " << stationarity << std::endl;
		std::cout << "constraints  eq: " << eq_con_values << std::endl;
		std::cout << "constraints ieq: " << ieq_con_values << std::endl;
	
		std::cout << "jac eq: " << nlp_->constraints_equality_jacobian(iterate) << std::endl;
	
		std::cout << "       dual  eq: " << iterate.dual_eq_values_ << std::endl;
		std::cout << "       dual ieq: " << iterate.dual_ieq_values_ << std::endl;
	}
	return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, iterate.dual_eq_values_, iterate.dual_ieq_values_);
}
double ResidualFunction::operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
	// first entry of rho(x,y,bary,mu):
	bool PRINT=true;
	std::vector<double> stationarity = nlp_->objective_gradient(iterate);
	
	std::shared_ptr<matrix_base_class> hessian = nlp_->lagrangian_hessian(iterate);
	
	std::vector<double> hessian_times_step = hessian->multiply(step.primal_values_);
	
	for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
		// mu*grad f + H*d:
		stationarity[stationarity_index] = iterate.penalty_parameter_*stationarity[stationarity_index] + hessian_times_step[stationarity_index];
		
	}
	
	std::shared_ptr<matrix_base_class> jacobian_eq = nlp_->constraints_equality_jacobian(iterate);
	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
	// std::cout << "about to do j eq" << std::endl;
	std::vector<double> linear_step_eq = jacobian_eq->multiply(step.primal_values_);
	
	// std::cout << "about to do j ieq" << std::endl;
	std::shared_ptr<matrix_base_class> jacobian_ieq = nlp_->constraints_inequality_jacobian(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
	std::vector<double> linear_step_ieq = jacobian_ieq->multiply(step.primal_values_);

	for (size_t eq_con_index=0; eq_con_index< nlp_->num_dual_eq(); ++eq_con_index) {
		con_values_eq[eq_con_index] += linear_step_eq[eq_con_index];
	}
	for (size_t ieq_con_index=0; ieq_con_index< nlp_->num_dual_ieq(); ++ieq_con_index) {
		con_values_ieq[ieq_con_index] += linear_step_ieq[ieq_con_index];
	}		
	return resid_helper(iterate, stationarity, con_values_eq, con_values_ieq, step.dual_eq_values_, step.dual_ieq_values_);
}