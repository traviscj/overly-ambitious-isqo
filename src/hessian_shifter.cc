
#include <cassert>

#include "step.hh"
// #include "iterate.hh"
#include "nlp.hh"
#include "nlp_state.hh"
#include "utilities.hh"
#include "subproblem.hh"
#include "solve_subproblem.hh"
#include "hessian_shifter.hh"


// HessianShifter::HessianShifter(Nlp &nlp) : FunctionWithNLPState(nlp), last_shift_(0.0) {
//     // TODO: support sparsity here.
//     // .... i think I could just have a HS::HS(SparseNlp &sparse_nlp) going.
//     // but the real question is: does it even matter?
//     // another question: maybe we should just accept a solver? (except currently they do not stick around.)
//     // or could also make the function not virtual void.
//     solve_qp_ = new SolveDenseQuadraticProgram(nlp);
// }

HessianShifter::HessianShifter(iSQOControlPanel &control, Nlp &nlp) : FunctionWithNLPState(control, nlp), solve_qp_(control, nlp), linear_model_reduction_(control, nlp), last_shift_(0.0), shifter_info("") {
    // solve_qp_ = new (nlp);
}

iSQOStep HessianShifter::operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem) {
    shifter_info.str("");
	bool PRINT=false;
	double shift_w0 = 1e-4;
	double shift_min = 1e-20;
	double shiftkwbarp = 100;
	double shiftkwp = 8;
	double shift_kwm = 1.0/3.0;
    //! \todo re-add shiftmax in hessian shifter.
    // double shiftmax = 1e40; 
	
	double current_shift = 0.0;
	
    shifter_info << "*** Entering the hessian shifter, penalty = " << iterate.get_penalty_parameter() << std::endl;
	// try to solve without regularization:
    // WHY would this line be necessary? (maybe so the matrix has the same sparsity structure?)
    subproblem.inc_regularization(current_shift, current_shift);
    int total_pivots=0;
	iSQOStep return_step = solve_qp_(subproblem, iterate);
    // std::cout << "hessianShifter :: operator() return_step (prior shift): " << return_step << std::endl;
	total_pivots += return_step.get_pivots();
	
	int current_regularization_steps = 0;
    // current_step_values has the length of the QP, because we want to multiply it with the QP hessian.
    // but we only copy the first (# NLP var), because we want it to be the reduction in the original space.
    // TODO check that with Andreas or Frank or Daniel.
	std::vector<double> current_step_values(subproblem.num_primal());
	for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
		current_step_values[primal_index] = return_step.get_primal_value(primal_index);
	}
	
	std::vector<double> hessian_times_step = subproblem.hessian_->multiply(current_step_values);
    if (PRINT) std::cout << "hessian step: " << subproblem.hessian_ << std::endl;
	double total = 0.0;
	for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
		total += hessian_times_step[primal_index] * return_step.get_primal_value(primal_index);
	}
	if (PRINT) std::cout << std::endl << "total: " << .5*total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
	if (PRINT) std::cout << "required: " << (.5*total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
    shifter_info << "**** Solve: "<< current_regularization_steps << "; Pivots: " << return_step.get_pivots() << "; Status: " << return_step.get_status() << " xxx: " << (.5*total < (1e-8*return_step.x_norm()*return_step.x_norm()))<< std::endl;
	// std::vector<double> hessian_step = subproblem.hessian_.multiply()
	while (        (return_step.get_status() != 0) 
                // || (return_step.x_norm() > 1e9) 
                || (.5*total < (1e-8*return_step.x_norm()*return_step.x_norm()))
                || linear_model_reduction_(iterate, return_step) < 0
                ) {
        double last_shift=current_shift;
		if (current_regularization_steps == 0) {
			if (last_shift_ == 0.0) {
				current_shift = shift_w0;
			} else {
				current_shift = fmax(shift_min, shift_kwm * last_shift_);
			}
		} else {
			if (last_shift_ == 0.0) {
				current_shift = shiftkwbarp * current_shift;
			} else {
				current_shift = shiftkwp*current_shift;
			}
		}
		subproblem.inc_regularization(current_shift, last_shift);
		return_step = solve_qp_(subproblem, iterate);
        int pivots_this_qp=return_step.get_pivots();
        // shifter_info << "**** Solve; Pivots: " << return_step.get_pivots() << "; Status: " << return_step.get_status() << std::endl;
        total_pivots += return_step.get_pivots();
        return_step.set_pivots(total_pivots);
		
		for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
			current_step_values[primal_index] = return_step.get_primal_value(primal_index);
		}
		if (PRINT) std::cout << "hessian: " << subproblem.hessian_ << std::endl;
		if (PRINT) std::cout << "return step: " << return_step.get_primal_values() << std::endl;
		hessian_times_step = subproblem.hessian_->multiply(current_step_values);
		if (PRINT) std::cout << "hessian step: " << hessian_times_step << std::endl;
		total = 0.0;
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
			total += hessian_times_step[primal_index] * return_step.get_primal_value(primal_index);
		}
		if (PRINT) std::cout << "current_shift: " << current_shift << std::endl;
		if (PRINT) std::cout << std::endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << std::endl;
		if (PRINT) std::cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << std::endl;
		
		++current_regularization_steps;

        shifter_info << "**** Solve: "<< current_regularization_steps << "; Pivots: " << pivots_this_qp << "; Status: " << return_step.get_status() << " xxx: " << (.5*total < (1e-8*return_step.x_norm()*return_step.x_norm()))<< std::endl;
		
		if (current_regularization_steps == 20){
			std::cout << "FAILURE!" << std::endl;
			assert(false);
			break;
		}
	}
    if (PRINT) std::cout << "SUCCESSFULLY SHIFTED!!" << current_regularization_steps << std::endl;
	last_shift_ = current_shift;
	return return_step;
}

double HessianShifter::get_last_shift() const {
	return last_shift_;
}
