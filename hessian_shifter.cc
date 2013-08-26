
#include <cassert>

#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_state.hh"
#include "utilities.hh"
#include "subproblem.hh"
#include "solve_subproblem.hh"
#include "hessian_shifter.hh"


HessianShifter::HessianShifter(Nlp &nlp) : FunctionWithNLPState(nlp), solve_qp_(nlp), last_shift_(0.0) {
	
}
iSQOStep HessianShifter::operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem) {
	bool PRINT=false;
	double shift_w0 = 1e-4;
	double shift_min = 1e-20;
	double shiftkwbarp = 100;
	double shiftkwp = 8;
	double shift_kwm = 1.0/3.0;
	double shiftmax = 1e40;
	
	double current_shift = 0.0;
	
	// try to solve without regularization:
	subproblem.inc_regularization(current_shift);
	iSQOStep return_step = solve_qp_(subproblem);
	
	
	int current_regularization_steps = 0;
	std::vector<double> current_step_values(nlp_->num_primal() + 2*nlp_->num_dual());
	for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
		current_step_values[primal_index] = return_step.primal_values_[primal_index];
	}
	
	// matrix hessian = nlp_->lagrangian_hessian(iterate);
	std::vector<double> hessian_step = subproblem.hessian_.multiply(current_step_values);
	if (PRINT) std::cout << "hessian step: " << subproblem.hessian_ << std::endl;
	double total = 0.0;
	for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
		total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
	}
	if (PRINT) std::cout << std::endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
	if (PRINT) std::cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
	// std::vector<double> hessian_step = subproblem.hessian_.multiply()
	while ((return_step.status_ != 0) || (return_step.x_norm() > 1e9) || (total < (1e-8*return_step.x_norm()*return_step.x_norm()))) {
		if (current_regularization_steps == 0) {
			if (last_shift_ == 0.0) {
				current_shift = shift_w0;
			} else {
				current_shift = std::max(shift_min, shift_kwm * last_shift_);
			}
		} else {
			if (last_shift_ == 0.0) {
				current_shift = shiftkwbarp * current_shift;
			} else {
				current_shift = shiftkwp*current_shift;
			}
		}
		subproblem.inc_regularization(current_shift);
		return_step = solve_qp_(subproblem);
		
		for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
			current_step_values[primal_index] = return_step.primal_values_[primal_index];
		}
		if (PRINT) std::cout << "hessian: " << subproblem.hessian_ << std::endl;
		if (PRINT) std::cout << "return step: " << return_step.primal_values_ << std::endl;
		hessian_step = subproblem.hessian_.multiply(current_step_values);
		if (PRINT) std::cout << "hessian step: " << hessian_step << std::endl;
		total = 0.0;
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
			total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
		}
		if (PRINT) std::cout << "current_shift: " << current_shift << std::endl;
		if (PRINT) std::cout << std::endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
		if (PRINT) std::cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
		
		++current_regularization_steps;
		
		if (current_regularization_steps == 20){
			std::cout << "FAILURE!" << std::endl;
			assert(false);
			break;
		}
	}
	last_shift_ = current_shift;
	return return_step;
}

double HessianShifter::get_last_shift() const {
	return last_shift_;
}
