
#include <iostream>
#include <vector>
#include <cmath>

#include "utilities.hh"
#include "constraint_violation.hh"


ConstraintViolationFunction::ConstraintViolationFunction(Nlp &nlp) : FunctionWithNLPState(nlp){
	// std::cout << "-- Initializing l1 violation function." << std::endl;
}
double ConstraintViolationFunction::operator()(const iSQOIterate &iterate) const {
    bool PRINT = false;
    if (PRINT) std::cout << "ConstraintViolationFunction @ iterate serial: " << iterate.get_serial() << std::endl;
	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
	
    if (PRINT) std::cout << "con_values_eq : " << con_values_eq << std::endl;
    if (PRINT) std::cout << "con_values_ieq: " << con_values_ieq << std::endl;
    
	return two_vectors(con_values_eq, con_values_ieq);
}
double ConstraintViolationFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
	bool PRINT = false;
    if (PRINT) std::cout << "ConstraintViolationFunction @ iterate serial: " << iterate.get_serial() << " + STEP: "<< std::endl;

	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);

	// J(x)*d
	if (PRINT) std::cout << "constraintviolationfunction e 0: " << con_values_eq << std::endl;
	std::shared_ptr<matrix_base_class> jac_ce = nlp_->constraints_equality_jacobian(iterate);
    if (PRINT) std::cout << "jac_ce: " << jac_ce << std::endl;
    if (PRINT) std::cout << "step: " << step.get_primal_values() << std::endl;
    std::vector<double> con_values_eq_mult_result = jac_ce->multiply(step.get_primal_values());
    for (size_t row=0; row < nlp_->num_dual_eq(); ++row) {
        con_values_eq[row] += con_values_eq_mult_result[row];
    }
	if (PRINT) std::cout << "constraintviolationfunction e 1: " << con_values_eq << std::endl;
	
	// \bar{J}(x)*d
	if (PRINT) std::cout << "constraintviolationfunction i 0: " << con_values_ieq << std::endl;
	std::shared_ptr<matrix_base_class> jac_ci = nlp_->constraints_inequality_jacobian(iterate);
    if (PRINT) std::cout << "jac_ci: " << jac_ci << std::endl;
    if (PRINT) std::cout << "step: " << step.get_primal_values() << std::endl;
    // std::vector<double> x(nlp_->num_primal());
    // x.assign(step.get_primal_values().begin(), step.get_primal_values().end());
    // std::vector<double> con_values_ieq_mult_result = jac_ci->multiply(x);
    std::vector<double> con_values_ieq_mult_result = jac_ci->multiply(step.get_primal_values());
    
    // std::cout << "con_values_ieq_mult_result     : " << con_values_ieq_mult_result << std::endl;
    // std::cout << "con_values_ieq_mult_result_orig: " << con_values_ieq_mult_result_orig << std::endl;
    if (PRINT) std::cout << "constraintviolationfunction step product: " << con_values_ieq_mult_result << std::endl;
    for (size_t row=0; row < nlp_->num_dual_ieq(); ++row) {
        con_values_ieq[row] = +con_values_ieq[row] + con_values_ieq_mult_result[row];
    }
	if (PRINT) std::cout << "constraintviolationfunction i 1: " << con_values_ieq << std::endl;
	
	return two_vectors(con_values_eq, con_values_ieq);
}

double ConstraintViolationFunction::two_vectors(const std::vector<double> &eq_items, const std::vector<double> &ieq_items) const {
    bool PRINT = false;
	double eq_violation = 0.0;
	for (size_t dual_eq_index=0; dual_eq_index<eq_items.size(); ++dual_eq_index) {
        if (PRINT) std::cout << "eq_violation pre step dual_eq_index=" << dual_eq_index << ": " << eq_violation ;
		eq_violation += fabs(eq_items[dual_eq_index]);
        if (PRINT) std::cout << "; eq_violation @ step dual_eq_index=" << dual_eq_index << ": " << eq_violation << std::endl;
	}
	double ieq_violation = 0.0;
	for (size_t dual_ieq_index=0; dual_ieq_index<ieq_items.size(); ++dual_ieq_index) {
		// if (ieq_items[i] > 0){
        if (PRINT) std::cout << "ieq_violation pre step dual_ieq_index=" << dual_ieq_index << ": " << ieq_violation << "; " << ieq_items[dual_ieq_index];
		ieq_violation += bracket_plus(ieq_items[dual_ieq_index]);
        if (PRINT) std::cout << "; ieq_violation @ step dual_ieq_index=" << dual_ieq_index << ": " << ieq_violation << std::endl;
		// }
	}
    if (PRINT) std::cout << "constraintviolationfunction  eq_violation: " << eq_violation << std::endl;
	if (PRINT) std::cout << "constraintviolationfunction ieq_violation: " << ieq_violation << std::endl;
    
	return eq_violation + ieq_violation;
}
