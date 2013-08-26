
#include <iostream>
#include <vector>

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

	if (PRINT) std::cout << "constraintviolationfunction e 0: " << con_values_eq << std::endl;
	// J(x)*d
	matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
	if (PRINT) std::cout << jac_ce;
	for (size_t row=0; row<iterate.num_dual_eq_; ++row) {
		for (size_t col=0; col<iterate.num_primal_; ++col) {
			con_values_eq[row] += jac_ce.get(row,col)*step.primal_values_[col];
		}
	}
	if (PRINT) std::cout << "constraintviolationfunction e 1: " << con_values_eq << std::endl;
	
	// \bar{J}(x)*d
	if (PRINT) std::cout << "constraintviolationfunction i 0: " << con_values_ieq << std::endl;
	matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
	// if (PRINT) jac_ci.print();
	for (size_t row=0; row<iterate.num_dual_ieq_; ++row) {
		for (size_t col=0; col<iterate.num_primal_; ++col) {
			con_values_ieq[row] += jac_ci.get(row,col)*step.primal_values_[col];
		}
	}
	if (PRINT) std::cout << "constraintviolationfunction i 1: " << con_values_ieq << std::endl;
	
	return two_vectors(con_values_eq, con_values_ieq);
}

double ConstraintViolationFunction::two_vectors(const std::vector<double> &eq_items, const std::vector<double> &ieq_items) const {
	double eq_violation = 0.0;
	for (size_t dual_eq_index=0; dual_eq_index<eq_items.size(); ++dual_eq_index) {
		eq_violation += abs(eq_items[dual_eq_index]);
	}
	double ieq_violation = 0.0;
	for (size_t dual_ieq_index=0; dual_ieq_index<ieq_items.size(); ++dual_ieq_index) {
		// if (ieq_items[i] > 0){
		ieq_violation += bracket_plus(ieq_items[dual_ieq_index]);
		// }
	}
	return eq_violation + ieq_violation;
}
