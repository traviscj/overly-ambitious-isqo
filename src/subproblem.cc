
#include <cmath>

#include "subproblem.hh"
#include "utilities.hh"

iSQOQuadraticSubproblem::iSQOQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate) :
			FunctionWithNLPState(nlp),
			num_qp_variables_(nlp.num_primal() + 2*nlp.num_dual()), 
			num_qp_constraints_(nlp.num_dual()),
			num_nlp_variables_(nlp.num_primal()), 
			num_nlp_constraints_eq_(nlp.num_dual_eq()), 
			num_nlp_constraints_ieq_(nlp.num_dual_ieq()),
			hessian_(num_qp_variables_, num_qp_variables_), 
			unshifted_hessian_diagonal_(num_qp_variables_), 
			first_shift_(false),
			jacobian_(num_qp_constraints_, num_qp_variables_), 
			gradient_(num_qp_variables_),
			lower_bound_(num_qp_variables_), 
			upper_bound_(num_qp_variables_), 
			jacobian_lower_bound_(num_qp_constraints_), 
			jacobian_upper_bound_(num_qp_constraints_),
			nlp_hessian_(num_nlp_variables_, num_nlp_variables_),
			nlp_eq_jacobian_(num_nlp_constraints_eq_, num_nlp_variables_),
			nlp_ieq_jacobian_(num_nlp_constraints_ieq_, num_nlp_variables_),
			nlp_objective_gradient_(num_nlp_variables_)
 {
	nlp_objective_gradient_ = nlp_->objective_gradient(iterate);
	
    // setup matrix data
    setup_matrix_data(iterate);
    
	// NLP gradient copied over
	for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index)
		gradient_[primal_index] = iterate.penalty_parameter_*nlp_objective_gradient_[primal_index];
	// equality slacks both have positive:
	for (size_t dual_eq_index=0; dual_eq_index<iterate.num_dual_eq_; ++dual_eq_index) {
		gradient_[iterate.num_primal_+dual_eq_index] = 1.0;
		gradient_[iterate.num_primal_+iterate.num_dual_eq_+dual_eq_index] = 1.0;
	}
	// inequality slacks only penalize the positive parts:
	for (size_t dual_ieq_index=0; dual_ieq_index<iterate.num_dual_ieq_; ++dual_ieq_index) {
		gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+dual_ieq_index] = 1.0;
		gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+dual_ieq_index] = 0.0;
	}
	
	std::vector<double> con_values_eq=nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);
	
	for (size_t eq_index=0; eq_index<iterate.num_dual_eq_; ++eq_index) {
		jacobian_lower_bound_[eq_index] = -con_values_eq[eq_index];
		jacobian_upper_bound_[eq_index] = -con_values_eq[eq_index];
	}
	for (size_t ieq_index=0; ieq_index<iterate.num_dual_ieq_; ++ieq_index) {
		jacobian_lower_bound_[iterate.num_dual_eq_+ieq_index] = -INFINITY;
		jacobian_upper_bound_[iterate.num_dual_eq_+ieq_index] = -con_values_ieq[ieq_index];
	}
	for (size_t variable_index=0; variable_index < iterate.num_primal_; ++variable_index) {
		lower_bound_[variable_index] = -INFINITY;
		upper_bound_[variable_index] = +INFINITY;
	}
	for (size_t variable_index=0; variable_index < 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_); ++variable_index) {
		lower_bound_[iterate.num_primal_+variable_index] = 0.0;
		upper_bound_[iterate.num_primal_+variable_index] = +INFINITY;
	}
}

void iSQOQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate) {
	nlp_hessian_ = nlp_->lagrangian_hessian(iterate);
	nlp_eq_jacobian_ = nlp_->constraints_equality_jacobian(iterate);
	nlp_ieq_jacobian_ = nlp_->constraints_inequality_jacobian(iterate);

	for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
		for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
			jacobian_.set(eq_constraint_index,variables, nlp_eq_jacobian_.get(eq_constraint_index,variables));
		}
		jacobian_.set(eq_constraint_index,nlp_->num_primal()+eq_constraint_index, -1.0);
		jacobian_.set(eq_constraint_index,nlp_->num_primal()+iterate.num_dual_eq_+iterate.num_dual_ieq_+eq_constraint_index, +1.0);
	}
			
	for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
		for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,variables, nlp_ieq_jacobian_.get(ieq_constraint_index,variables));
		}
		jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+iterate.num_dual_eq_+ieq_constraint_index, -1.0);
		jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+iterate.num_dual_eq_+iterate.num_dual_ieq_+iterate.num_dual_eq_+ieq_constraint_index, +1.0);
	}
    
	for (size_t r=0; r<iterate.num_primal_; ++r) {
		for (size_t c=0; c<iterate.num_primal_; ++c) {
			hessian_.set(r,c, nlp_hessian_.get(r,c));
		}
	}
}



void iSQOQuadraticSubproblem::inc_regularization(double hessian_shift) {
	bool PRINT = false;
	// make a copy of the diagonal entries, if we haven't already got one...
	if (! first_shift_) {
		if (PRINT) std::cout << "making the first-time-only copy..." << std::endl;
		for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
				unshifted_hessian_diagonal_[diagonal_entry] = hessian_.get(diagonal_entry,diagonal_entry);
		}
		first_shift_ = true;
	}
	
	for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
			hessian_.set(diagonal_entry,diagonal_entry, unshifted_hessian_diagonal_[diagonal_entry] + hessian_shift);
	}
	for (size_t diagonal_entry=num_nlp_variables_; diagonal_entry < num_qp_variables_; ++diagonal_entry) {
		hessian_.set(diagonal_entry, diagonal_entry, hessian_shift);
	}
}

void iSQOQuadraticSubproblem::print() const {
	std::cout << std::endl;
	std::cout << "H = [";
	std::cout << hessian_;
	std::cout << "];" << std::endl;
	std::cout << "A = [";
	std::cout << jacobian_;
	std::cout << "];" << std::endl;
	
	std::cout << "g = " << gradient_ << "';" << std::endl;
	std::cout << "lb = " << lower_bound_ << "';" << std::endl;
	std::cout << "ub = " << upper_bound_ << "';" << std::endl;
	std::cout << "lbA = " << jacobian_lower_bound_ << "';" << std::endl;
	std::cout << "ubA = " << jacobian_upper_bound_ << "';" << std::endl;
}

iSQOSparseQuadraticSubproblem::iSQOSparseQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate) :
    iSQOQuadraticSubproblem(nlp, iterate),
    hessian_sparse_(nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(), 0),
    jacobian_sparse_(nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(), 0),
    nlp_hessian_sparse_(nlp.num_primal() + 2*nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(),500),
    nlp_eq_jacobian_sparse_(nlp.num_dual_eq(), nlp.num_primal() + 2*nlp.num_dual(),500),
    nlp_ieq_jacobian_sparse_(nlp.num_dual_ieq(), nlp.num_primal() + 2*nlp.num_dual(),500)
    {
        
        setup_matrix_data_sparse(iterate);
}

void iSQOSparseQuadraticSubproblem::setup_matrix_data_sparse(const iSQOIterate &iterate) {
	
    // std::cout << "iSQOQuadraticSubproblem::setup_matrix_data" << std::endl;
    
    // JACOBIAN PART
	nlp_eq_jacobian_sparse_ = nlp_->constraints_equality_jacobian_sparse(iterate);
    // std::cout << "nlp_eq_jacobian_sparse_: " << nlp_eq_jacobian_sparse_ << std::endl;
	nlp_ieq_jacobian_sparse_ = nlp_->constraints_inequality_jacobian_sparse(iterate);
    // std::cout << "nlp_ieq_jacobian_sparse_: " << nlp_ieq_jacobian_sparse_ << std::endl;
	
    sparse_matrix jacobian_eq_ieq = vertical(nlp_eq_jacobian_sparse_, nlp_ieq_jacobian_sparse_);
    // std::cout << "jacobian_eq_ieq: " << jacobian_eq_ieq << std::endl;
    
    sparse_matrix jacobian_p_slack(iterate.num_dual_eq_+iterate.num_dual_ieq_, -1.0);
    // std::cout << "jacobian_p_slack: " << jacobian_p_slack << std::endl;
    sparse_matrix jacobian_n_slack(iterate.num_dual_eq_+iterate.num_dual_ieq_, +1.0);
    // std::cout << "jacobian_n_slack: " << jacobian_n_slack << std::endl;
    sparse_matrix jacobian_slacks = horizontal(jacobian_p_slack, jacobian_n_slack);
    // std::cout << "jacobian_slacks: " << jacobian_slacks << std::endl;
    
    jacobian_sparse_ = horizontal(jacobian_eq_ieq,jacobian_slacks);
    // std::cout << "jacobian_sparse: " << jacobian_sparse_ << std::endl;
    
    // LAGRANGIAN PART
    nlp_hessian_sparse_ = nlp_->lagrangian_hessian_sparse(iterate);
    sparse_matrix hess_left_bottom(2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), num_nlp_variables_, 0);    // null matrix for lower left
    sparse_matrix hess_right_top(num_nlp_variables_, 2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), 0);  // null matrix for upper right.
    sparse_matrix hess_right_bottom(2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), 0.0);            // 0.0*I_{2*num_qp_con} for lower right. (fixes bug in qpOASES...)
    // std::cout << "hess_right_bottom: " << hess_right_bottom << std::endl;
    
    sparse_matrix hess_left = vertical(nlp_hessian_sparse_, hess_left_bottom);
    // std::cout << "hess_left: " << hess_left << std::endl;
    sparse_matrix hess_right = vertical(hess_right_top, hess_right_bottom);
    // std::cout << "hess_right: " << hess_right << std::endl;
    hessian_sparse_ = horizontal(hess_left, hess_right);
    // std::cout << "hessian_sparse_: " << hessian_sparse_ << std::endl;
    
}
void iSQOSparseQuadraticSubproblem::inc_regularization(double hessian_shift) {
	bool PRINT = false;
    // make a copy of the diagonal entries, if we haven't already got one...
    if (! first_shift_) {
        if (PRINT) std::cout << "making the first-time-only copy..." << std::endl;
        for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
                unshifted_hessian_diagonal_[diagonal_entry] = hessian_.get(diagonal_entry,diagonal_entry);
        }
        first_shift_ = true;
    }
    hessian_shift_ = hessian_shift;
    for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
            hessian_.set(diagonal_entry,diagonal_entry, unshifted_hessian_diagonal_[diagonal_entry] + hessian_shift);
    }
    for (size_t diagonal_entry=num_nlp_variables_; diagonal_entry < num_qp_variables_; ++diagonal_entry) {
        hessian_.set(diagonal_entry, diagonal_entry, hessian_shift);
    }
}
