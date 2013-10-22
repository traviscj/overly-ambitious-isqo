// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#include <iostream>
#include <vector>
#include <cmath>

#include "utilities.hh"
#include "constraint_violation.hh"


ConstraintViolationFunction::ConstraintViolationFunction(Nlp &nlp) : FunctionWithNLPState(nlp){
	// std::cout << "-- Initializing l1 violation function." << std::endl;
}
double ConstraintViolationFunction::operator()(const iSQOIterate &iterate) const {
	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
	
	return two_vectors(con_values_eq, con_values_ieq);
}
double ConstraintViolationFunction::operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
	bool PRINT = false;

	std::vector<double> con_values_eq = nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);

	// J(x)*d
    // if (PRINT) std::cout << "jac_ce: " << jac_ce << std::endl;
    // if (PRINT) std::cout << "step: " << step.get_primal_values() << std::endl;
	if (PRINT) std::cout << "constraintviolationfunction e 0: " << con_values_eq << std::endl;
	std::shared_ptr<matrix_base_class> jac_ce = nlp_->constraints_equality_jacobian(iterate);
    std::vector<double> con_values_eq_mult_result = jac_ce->multiply(step.get_primal_values());
    for (size_t row=0; row < nlp_->num_dual_eq(); ++row) {
        con_values_eq[row] += con_values_eq_mult_result[row];
    }
	if (PRINT) std::cout << "constraintviolationfunction e 1: " << con_values_eq << std::endl;
	
	// \bar{J}(x)*d
    // if (PRINT) std::cout << "jac_ci: " << jac_ci << std::endl;
    // if (PRINT) std::cout << "step: " << step.get_primal_values() << std::endl;
	if (PRINT) std::cout << "constraintviolationfunction i 0: " << con_values_ieq << std::endl;
	std::shared_ptr<matrix_base_class> jac_ci = nlp_->constraints_inequality_jacobian(iterate);
    std::vector<double> con_values_ieq_mult_result = jac_ci->multiply(step.get_primal_values());
    for (size_t row=0; row < nlp_->num_dual_ieq(); ++row) {
        con_values_ieq[row] += con_values_ieq_mult_result[row];
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
		ieq_violation += bracket_plus(ieq_items[dual_ieq_index]);
		// }
	}
    if (PRINT) std::cout << "constraintviolationfunction  eq_violation: " << eq_violation << std::endl;
	if (PRINT) std::cout << "constraintviolationfunction ieq_violation: " << ieq_violation << std::endl;
    
	return eq_violation + ieq_violation;
}
