
#include <cassert>
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
			hessian_(NULL), 
			unshifted_hessian_diagonal_(num_qp_variables_), 
			first_shift_(false),
			jacobian_(NULL), 
			gradient_(num_qp_variables_),
			lower_bound_(num_qp_variables_), 
			upper_bound_(num_qp_variables_), 
			jacobian_lower_bound_(num_qp_constraints_), 
			jacobian_upper_bound_(num_qp_constraints_),
			nlp_hessian_(NULL),
			nlp_eq_jacobian_(NULL),
			nlp_ieq_jacobian_(NULL),
			nlp_objective_gradient_(num_nlp_variables_),
            iterate_pointer(&iterate)
 {
	nlp_objective_gradient_ = nlp_->objective_gradient(iterate);
	    
    // setup matrix data
    setup_matrix_data(iterate, nlp_->constraints_equality_jacobian(iterate), nlp_->constraints_inequality_jacobian(iterate), nlp_->lagrangian_hessian(iterate));
    
	// NLP gradient copied over
	for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
		gradient_[primal_index] = iterate.get_penalty_parameter()*nlp_objective_gradient_[primal_index];
	// equality slacks both have positive:
	for (size_t dual_eq_index=0; dual_eq_index<nlp_->num_dual_eq(); ++dual_eq_index) {
		gradient_[nlp_->num_primal()+dual_eq_index] = 1.0;
		gradient_[nlp_->num_primal()+nlp_->num_dual_eq()+dual_eq_index] = 1.0;
	}
	// inequality slacks only penalize the positive parts:
	for (size_t dual_ieq_index=0; dual_ieq_index<nlp_->num_dual_ieq(); ++dual_ieq_index) {
		gradient_[nlp_->num_primal()+2*nlp_->num_dual_eq()+dual_ieq_index] = 1.0;
		gradient_[nlp_->num_primal()+2*nlp_->num_dual_eq()+nlp_->num_dual_ieq()+dual_ieq_index] = 0.0;
	}
	
	std::vector<double> con_values_eq=nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);

    // formulation 1: group p - X^+
    size_t BOUND_REFORMULATION_MODE = 2;
    if (BOUND_REFORMULATION_MODE == 1) {
        assert(false); // this formulation is WRONG!
        for (size_t eq_index=0; eq_index<nlp_->num_dual_eq(); ++eq_index) {
            jacobian_lower_bound_[eq_index] = 0.0; // was: -con_values_eq[eq_index];
            jacobian_upper_bound_[eq_index] = 0.0; // was: -con_values_eq[eq_index];
            lower_bound_[nlp_->num_primal() + eq_index] = -bracket_plus(con_values_eq[eq_index]);
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + eq_index] = -bracket_minus(con_values_eq[eq_index]);
            upper_bound_[nlp_->num_primal() + eq_index] = INFINITY;
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + eq_index] = INFINITY;
        }
        for (size_t ieq_index=0; ieq_index<nlp_->num_dual_ieq(); ++ieq_index) {
            jacobian_lower_bound_[nlp_->num_dual_eq()+ieq_index] = 0.0; // was: -INFINITY;
            jacobian_upper_bound_[nlp_->num_dual_eq()+ieq_index] = 0.0; // was: -con_values_ieq[ieq_index];
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + ieq_index] = -bracket_plus(con_values_ieq[ieq_index]);
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + nlp_->num_dual_eq() + ieq_index] = -bracket_minus(con_values_ieq[ieq_index]);
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + ieq_index] = INFINITY;
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + nlp_->num_dual_eq() + ieq_index] = INFINITY;            
        }
    } else if (BOUND_REFORMULATION_MODE == 2) { 
        // formulation two: grouping p + X^+:
    	for (size_t eq_index=0; eq_index<nlp_->num_dual_eq(); ++eq_index) {
    		jacobian_lower_bound_[eq_index] = 0.0; // was: -con_values_eq[eq_index];
    		jacobian_upper_bound_[eq_index] = 0.0; // was: -con_values_eq[eq_index];
            lower_bound_[nlp_->num_primal() + eq_index] = bracket_minus(con_values_eq[eq_index]);
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + eq_index] = bracket_plus(con_values_eq[eq_index]);
            upper_bound_[nlp_->num_primal() + eq_index] = INFINITY;
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + eq_index] = INFINITY;
    	}
    	for (size_t ieq_index=0; ieq_index<nlp_->num_dual_ieq(); ++ieq_index) {
    		jacobian_lower_bound_[nlp_->num_dual_eq()+ieq_index] = 0.0; // was: -INFINITY;
    		jacobian_upper_bound_[nlp_->num_dual_eq()+ieq_index] = 0.0; // was: -con_values_ieq[ieq_index];
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + ieq_index] = bracket_minus(con_values_ieq[ieq_index]);
            lower_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + nlp_->num_dual_eq() + ieq_index] = bracket_plus(con_values_ieq[ieq_index]);
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + ieq_index] = INFINITY;
            upper_bound_[nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq() + nlp_->num_dual_eq() + ieq_index] = INFINITY;
        
    	}
    } else {
        for (size_t eq_index=0; eq_index<nlp_->num_dual_eq(); ++eq_index) {
            jacobian_lower_bound_[eq_index] = -con_values_eq[eq_index];
            jacobian_upper_bound_[eq_index] = -con_values_eq[eq_index];
        }
        for (size_t ieq_index=0; ieq_index<nlp_->num_dual_ieq(); ++ieq_index) {
            jacobian_lower_bound_[nlp_->num_dual_eq()+ieq_index] = -INFINITY;
            jacobian_upper_bound_[nlp_->num_dual_eq()+ieq_index] = -con_values_ieq[ieq_index];
        }
        for (size_t variable_index=0; variable_index < 2*(nlp_->num_dual_eq() + nlp_->num_dual_ieq()); ++variable_index) {
            lower_bound_[nlp_->num_primal()+variable_index] = 0.0;
            upper_bound_[nlp_->num_primal()+variable_index] = +INFINITY;
        }
    }
	for (size_t variable_index=0; variable_index < nlp_->num_primal(); ++variable_index) {
		lower_bound_[variable_index] = -INFINITY;
		upper_bound_[variable_index] = +INFINITY;
	}
	for (size_t variable_index=0; variable_index < 2*(nlp_->num_dual_eq() + nlp_->num_dual_ieq()); ++variable_index) {
        // lower_bound_[nlp_->num_primal()+variable_index] = 0.0;
        // upper_bound_[nlp_->num_primal()+variable_index] = +INFINITY;
	}
    // std::cout << "jacobian_lower_bound_: " << jacobian_lower_bound_ << std::endl;
    // std::cout << "jacobian_upper_bound_: " << jacobian_upper_bound_ << std::endl;
    // std::cout << "lower_bound_: " << lower_bound_ << std::endl;
    // std::cout << "upper_bound_: " << upper_bound_ << std::endl;
    
}

void iSQOQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<dense_matrix> nlp_eq_jacobian, std::shared_ptr<dense_matrix> nlp_ieq_jacobian, std::shared_ptr<dense_matrix> nlp_hessian) {
	
    nlp_eq_jacobian_ = nlp_eq_jacobian;
    nlp_ieq_jacobian_ = nlp_ieq_jacobian;
    nlp_hessian_ = nlp_hessian;
    
    // jacobian_ = std::shared_ptr<dense_matrix>(new dense_matrix(num_qp_constraints_, num_qp_variables_));
    std::shared_ptr<dense_matrix> jacobian_tmp(new dense_matrix(num_qp_constraints_, num_qp_variables_));
    
	for (size_t eq_constraint_index=0; eq_constraint_index < nlp_->num_dual_eq(); ++eq_constraint_index) {
		for (size_t variables=0; variables<nlp_->num_primal(); ++variables) {
			jacobian_tmp->set(eq_constraint_index,variables, nlp_eq_jacobian->get(eq_constraint_index,variables));
		}
		jacobian_tmp->set(eq_constraint_index,nlp_->num_primal()+eq_constraint_index, -1.0);
		jacobian_tmp->set(eq_constraint_index,nlp_->num_primal()+nlp_->num_dual_eq()+nlp_->num_dual_ieq()+eq_constraint_index, +1.0);
	}
			
	for (size_t ieq_constraint_index=0; ieq_constraint_index < nlp_->num_dual_ieq(); ++ieq_constraint_index) {
		for (size_t variables=0; variables<nlp_->num_primal(); ++variables) {
			jacobian_tmp->set(nlp_->num_dual_eq()+ieq_constraint_index,variables, nlp_ieq_jacobian->get(ieq_constraint_index,variables));
		}
		jacobian_tmp->set(nlp_->num_dual_eq()+ieq_constraint_index, nlp_->num_primal()+nlp_->num_dual_eq()+ieq_constraint_index, -1.0);
		jacobian_tmp->set(nlp_->num_dual_eq()+ieq_constraint_index, nlp_->num_primal()+nlp_->num_dual_eq()+nlp_->num_dual_ieq()+nlp_->num_dual_eq()+ieq_constraint_index, +1.0);
	}
    jacobian_ = jacobian_tmp;
    
    // hessian_ = std::shared_ptr<dense_matrix>(new dense_matrix(num_qp_variables_, num_qp_variables_));
    std::shared_ptr<dense_matrix> hessian_tmp(new dense_matrix(num_qp_variables_, num_qp_variables_));
	for (size_t r=0; r<nlp_->num_primal(); ++r) {
		for (size_t c=0; c<nlp_->num_primal(); ++c) {
			hessian_tmp->set(r,c, nlp_hessian->get(r,c));
		}
	}
    hessian_ = hessian_tmp;
}



void iSQOQuadraticSubproblem::inc_regularization(double hessian_shift, double last_shift) {
	bool PRINT = false;
    // std::cout << "hessian_shift: " << hessian_shift << "; last_shift: " << last_shift << std::endl;
    hessian_->regularize(hessian_shift, last_shift);
    // std::cout << "\n\n\n nlp_hessian_ pre: " << nlp_hessian_ << std::endl;
    nlp_hessian_->regularize(hessian_shift, last_shift);
    // std::cout << "\n\n\n nlp_hessian_ post: " << nlp_hessian_ << std::endl;
    
    // std::cout << "\n\n\n nlp hessian re-evaluation :" << nlp_->lagrangian_hessian(*iterate_pointer) << std::endl;
    // std::shared_ptr< identity(NUM_FUCKING_VARIABLES, (hessian_shift - last_shift));
    // hessian_ = sum(hessian_, )
}

std::ostream &iSQOQuadraticSubproblem::print(std::ostream &os) const {
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

// iSQOSparseQuadraticSubproblem::iSQOSparseQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate) :
//     iSQOQuadraticSubproblem(nlp, iterate),
//     hessian_sparse_(new sparse_matrix(nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(), 0)),
//     jacobian_sparse_(new sparse_matrix(nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(), 0)),
//     nlp_hessian_sparse_(new sparse_matrix(nlp.num_primal() + 2*nlp.num_dual(), nlp.num_primal() + 2*nlp.num_dual(),500)),
//     nlp_eq_jacobian_sparse_(new sparse_matrix(nlp.num_dual_eq(), nlp.num_primal() + 2*nlp.num_dual(),500)),
//     nlp_ieq_jacobian_sparse_(new sparse_matrix(nlp.num_dual_ieq(), nlp.num_primal() + 2*nlp.num_dual(),500))
//     {
//         
//         setup_matrix_data_sparse(iterate);
// }

void iSQOQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<matrix_base_class> nlp_eq_jacobian, std::shared_ptr<matrix_base_class> nlp_ieq_jacobian, std::shared_ptr<matrix_base_class> nlp_hessian){
    std::shared_ptr<dense_matrix> attempt_dense = std::dynamic_pointer_cast< dense_matrix  >(nlp_eq_jacobian);
    std::shared_ptr<sparse_matrix> attempt_sparse = std::dynamic_pointer_cast< sparse_matrix >(nlp_eq_jacobian);
    
    // std::cout << "setting up matrix data..." << std::endl;
    
    if (attempt_dense != NULL) {
        return setup_matrix_data(iterate, attempt_dense, std::dynamic_pointer_cast< dense_matrix  >(nlp_ieq_jacobian), std::dynamic_pointer_cast< dense_matrix  >(nlp_hessian));
    } else if (attempt_sparse != NULL) {
        return setup_matrix_data(iterate, attempt_sparse, std::dynamic_pointer_cast< sparse_matrix >(nlp_ieq_jacobian), std::dynamic_pointer_cast< sparse_matrix >(nlp_hessian));
    } else {
        std::cout << "wtf mate" << std::endl;
        throw("wtf mate");
    }
}


void iSQOQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<sparse_matrix> nlp_eq_jacobian, std::shared_ptr<sparse_matrix> nlp_ieq_jacobian, std::shared_ptr<sparse_matrix> nlp_hessian) {
    // void iSQOSparseQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate) {
	
    // std::cout << "iSQOQuadraticSubproblem::setup_matrix_data" << std::endl;
	nlp_eq_jacobian_ = nlp_eq_jacobian;
    // std::cout << "nlp_eq_jacobian_sparse_: " << nlp_eq_jacobian_sparse_ << std::endl;
	nlp_ieq_jacobian_ = nlp_ieq_jacobian;
    // std::cout << "nlp_ieq_jacobian_sparse_: " << nlp_ieq_jacobian_sparse_ << std::endl;
	nlp_hessian_ = nlp_hessian;
    
    // JACOBIAN PART
    //     std::map<int, std::shared_ptr<matrix_base_class> >::iterator jacobian_it = prior_constructed_jacobians_.find(iterate.get_serial());
    //     if (jacobian_it != prior_constructed_jacobians_.end()) {jacobian_ = jacobian_it->second; std::cout << "JAC HIT\n";}
    // else {
        std::shared_ptr<sparse_matrix> jacobian_eq_ieq = nlp_eq_jacobian->vertical(nlp_ieq_jacobian);
        // std::cout << "jacobian_eq_ieq: " << jacobian_eq_ieq << std::endl;
    
        std::shared_ptr<sparse_matrix> jacobian_p_slack(new sparse_matrix(nlp_->num_dual_eq()+nlp_->num_dual_ieq(), -1.0));
        // std::cout << "jacobian_p_slack: " << jacobian_p_slack << std::endl;
        std::shared_ptr<sparse_matrix> jacobian_n_slack(new sparse_matrix(nlp_->num_dual_eq()+nlp_->num_dual_ieq(), +1.0));
        // std::cout << "jacobian_n_slack: " << jacobian_n_slack << std::endl;
        std::shared_ptr<sparse_matrix> jacobian_slacks = jacobian_p_slack->horizontal(jacobian_n_slack);
        // std::cout << "jacobian_slacks: " << jacobian_slacks << std::endl;
    
        jacobian_ = jacobian_eq_ieq->horizontal(jacobian_slacks);
        // prior_constructed_jacobians_[iterate.get_serial()] = jacobian_;
        // std::cout << "jacobian_sparse: " << jacobian_sparse_ << std::endl;
    // }
    
    // LAGRANGIAN PART
    //     std::map<int, std::shared_ptr<matrix_base_class> >::iterator hessian_it = prior_constructed_hessians_.find(iterate.get_serial());
    //     if (hessian_it != prior_constructed_hessians_.end()) {hessian_ = hessian_it->second; std::cout << "HESS HIT\n";}
    // else {
        std::shared_ptr<sparse_matrix> hess_left_bottom(new sparse_matrix(2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), num_nlp_variables_, 0));    // null matrix for lower left
        std::shared_ptr<sparse_matrix> hess_right_top(new sparse_matrix(num_nlp_variables_, 2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), 0));  // null matrix for upper right.
        std::shared_ptr<sparse_matrix> hess_right_bottom(new sparse_matrix(2*(num_nlp_constraints_eq_ + num_nlp_constraints_ieq_), 0.0));            // 0.0*I_{2*num_qp_con} for lower right. (fixes bug in qpOASES...)
        // std::cout << "hess_right_bottom: " << hess_right_bottom << std::endl;
    
        std::shared_ptr<sparse_matrix> hess_left = nlp_hessian->vertical(hess_left_bottom);
        // std::cout << "hess_left: " << hess_left << std::endl;
        std::shared_ptr<sparse_matrix> hess_right = hess_right_top->vertical(hess_right_bottom);
        // std::cout << "hess_right: " << hess_right << std::endl;
        hessian_ = hess_left->horizontal(hess_right);
        // prior_constructed_hessians_[iterate.get_serial()] = hessian_;
    // }
    // std::cout << "hessian_sparse_: " << hessian_sparse_ << std::endl;
    
}
// void iSQOSparseQuadraticSubproblem::inc_regularization(double hessian_shift) {
//     bool PRINT = false;
//     // TODO THIS IS REALLY REALLY REALLY BROKEN
//     // make a copy of the diagonal entries, if we haven't already got one...
//     // if (! first_shift_) {
//     //     if (PRINT) std::cout << "making the first-time-only copy..." << std::endl;
//     //     for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
//     //             unshifted_hessian_diagonal_[diagonal_entry] = hessian_->get(diagonal_entry,diagonal_entry);
//     //     }
//     //     first_shift_ = true;
//     // }
//     // hessian_shift_ = hessian_shift;
//     // for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
//     //         hessian_->set(diagonal_entry,diagonal_entry, unshifted_hessian_diagonal_[diagonal_entry] + hessian_shift);
//     // }
//     // for (size_t diagonal_entry=num_nlp_variables_; diagonal_entry < num_qp_variables_; ++diagonal_entry) {
//     //     hessian_->set(diagonal_entry, diagonal_entry, hessian_shift);
//     // }
// }
