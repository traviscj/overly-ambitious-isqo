
#include <cassert>
#include <cmath>

#include "subproblem.hh"
#include "utilities.hh"

iSQOQuadraticSubproblem::iSQOQuadraticSubproblem(iSQOControlPanel &control, Nlp &nlp, const iSQOIterate &iterate) :
			FunctionWithNLPState(control, nlp),
			num_qp_variables_(nlp.num_primal() + 2*nlp.num_dual_eq() + nlp.num_dual_ieq()), 
			num_qp_constraints_(nlp.num_dual()),
			num_nlp_variables_(nlp.num_primal()), 
			num_nlp_constraints_eq_(nlp.num_dual_eq()), 
			num_nlp_constraints_ieq_(nlp.num_dual_ieq()),
			hessian_(NULL), 
			unshifted_hessian_diagonal_(num_qp_variables_), 
            // first_shift_(false),
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

    // setup indices into these arrays...
    std::size_t start_eq_p = nlp_->num_primal();
    std::size_t start_eq_n = nlp_->num_primal() + nlp_->num_dual_eq() + nlp_->num_dual_ieq();
    std::size_t start_ieq_slack = nlp_->num_primal() + nlp_->num_dual_eq();
    	
	std::vector<double> con_values_eq=nlp_->constraints_equality(iterate);
	std::vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);

	for (size_t primal_index=0; primal_index < nlp_->num_primal(); ++primal_index) {
        // setup NLP gradient:
        gradient_[primal_index] = iterate.get_penalty_parameter()*nlp_objective_gradient_[primal_index];
        
        // we ignore variable bounds, since we've already converted them to penalized constraints.
		lower_bound_[primal_index] = -INFINITY;
		upper_bound_[primal_index] = +INFINITY;
	}
     for (size_t dual_eq_index=0; dual_eq_index<nlp_->num_dual_eq(); ++dual_eq_index) {
         // S\ell_1QP gradient with penalty
 		gradient_[start_eq_p + dual_eq_index] = 1.0;
 		gradient_[start_eq_n + dual_eq_index] = 1.0;
         
        // Equality constraints rolled into Jacobian bounds:
        jacobian_lower_bound_[dual_eq_index] = -con_values_eq[dual_eq_index];
        jacobian_upper_bound_[dual_eq_index] = -con_values_eq[dual_eq_index];
        
        // setup p bounds:
        lower_bound_[start_eq_p + dual_eq_index] = 0.0;
        upper_bound_[start_eq_p + dual_eq_index] = INFINITY;
        
        // setup n bounds:
        lower_bound_[start_eq_n + dual_eq_index] = 0.0;
        upper_bound_[start_eq_n + dual_eq_index] = INFINITY;
    }
    for (size_t dual_ieq_index=0; dual_ieq_index<nlp_->num_dual_ieq(); ++dual_ieq_index) {
        // only need to set up the relaxation variable, not the slack.
        
        // S\ell_1QP gradient for 
		gradient_[start_ieq_slack + dual_ieq_index] = 1.0;
        
        // Inequality constraints rolled into bounds:
        jacobian_lower_bound_[nlp_->num_dual_eq()+dual_ieq_index] = -INFINITY;
        jacobian_upper_bound_[nlp_->num_dual_eq()+dual_ieq_index] = -con_values_ieq[dual_ieq_index];
        
        // setup \bar{p} bounds:
        lower_bound_[start_ieq_slack + dual_ieq_index] = 0.0;
        upper_bound_[start_ieq_slack + dual_ieq_index] = INFINITY;
    }
    
	
    // std::cout << "jacobian_lower_bound_: " << jacobian_lower_bound_ << std::endl;
    // std::cout << "jacobian_upper_bound_: " << jacobian_upper_bound_ << std::endl;
    // std::cout << "lower_bound_: " << lower_bound_ << std::endl;
    // std::cout << "upper_bound_: " << upper_bound_ << std::endl;
    
}

void iSQOQuadraticSubproblem::setup_matrix_data(const iSQOIterate &iterate, std::shared_ptr<dense_matrix> nlp_eq_jacobian, std::shared_ptr<dense_matrix> nlp_ieq_jacobian, std::shared_ptr<dense_matrix> nlp_hessian) {
	assert(false); // didn't update this for resizes...
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
    // bool con_values_ieqPRINT = false;
    std::cout << "hessian_shift: " << hessian_shift << "; last_shift: " << last_shift << std::endl;
    hessian_->regularize(hessian_shift, last_shift);
    // std::cout << "\n\n\n nlp_hessian_ pre: " << nlp_hessian_ << std::endl;
    nlp_hessian_->regularize(hessian_shift, last_shift);
    // std::cout << "\n\n\n nlp_hessian_ post: " << nlp_hessian_ << std::endl;
    
    // std::cout << "\n\n\n nlp hessian re-evaluation :" << nlp_->lagrangian_hessian(*iterate_pointer) << std::endl;
    // std::shared_ptr< identity(NUM_FUCKING_VARIABLES, (hessian_shift - last_shift));
    // hessian_ = sum(hessian_, )
}

std::ostream &iSQOQuadraticSubproblem::print(std::ostream &os) const {
	os << std::endl;
	os << "mat = sparse(" << num_qp_variables_ << ", " << num_qp_variables_ << ");" << std::endl;
	os << hessian_;
    os << "H=mat;" << std::endl;
    os << "mat = sparse(" << num_qp_constraints_ << ", " << num_qp_variables_ << ");" << std::endl;
	os << jacobian_;
    os << "A=mat;" << std::endl;
    
    // std::cout << "];" << std::endl;
	
	os << "g = " << gradient_ << "';" << std::endl;
	os << "lb = " << lower_bound_ << "';" << std::endl;
	os << "ub = " << upper_bound_ << "';" << std::endl;
	os << "lbA = " << jacobian_lower_bound_ << "';" << std::endl;
	os << "ubA = " << jacobian_upper_bound_ << "';" << std::endl;
    return os;
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
    bool PRINT = false;
	
    // std::cout << "iSQOQuadraticSubproblem::setup_matrix_data" << std::endl;
	nlp_eq_jacobian_ = nlp_eq_jacobian;
    if (PRINT) std::cout << "nlp_eq_jacobian_: " << nlp_eq_jacobian_ << std::endl;
	nlp_ieq_jacobian_ = nlp_ieq_jacobian;
    if (PRINT) std::cout << "nlp_ieq_jacobian_: " << nlp_ieq_jacobian_ << std::endl;
	nlp_hessian_ = nlp_hessian;
    
    int n_eq = nlp_->num_dual_eq();
    int n_ieq = nlp_->num_dual_ieq();
    
    // JACOBIAN PART
        std::map<int, std::shared_ptr<matrix_base_class> >::iterator jacobian_it = prior_constructed_jacobians_.find(iterate.get_serial());
        if (jacobian_it != prior_constructed_jacobians_.end()) {jacobian_ = jacobian_it->second; std::cout << "JAC HIT\n";}
    else {
        std::shared_ptr<sparse_matrix> jacobian_eq_ieq = nlp_eq_jacobian->vertical(nlp_ieq_jacobian);
        if (PRINT) std::cout << "jacobian_eq_ieq: " << *jacobian_eq_ieq << std::endl;
    
        // Here, we're going to construct the slacks matrix:
        // / -I  0 +I \
        // \ 0  -I  0 /
        // where the null/identity matrices are referred to below by jacobian_slack_(block row)(block col).
        // Defining
        // matrices, by block row              block row--+     +--block column     #rows, #cols/scalar, num_nonzeros)
        // ======================                         v     v                   -----  -----         -----
        std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol0(sparse_matrix::diagonal_matrix( n_eq, n_eq, -1.0));
        std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol1(new sparse_matrix( n_eq, n_ieq, 0));
        std::shared_ptr<sparse_matrix> jacobian_slack_brow0_bcol2(sparse_matrix::diagonal_matrix( n_eq, n_eq,  +1.0) );
    
        std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol0(new sparse_matrix(n_ieq,  n_eq,      0));
        std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol1(sparse_matrix::diagonal_matrix(n_ieq, n_ieq, -1.0));
        std::shared_ptr<sparse_matrix> jacobian_slack_brow1_bcol2(new sparse_matrix(n_ieq,  n_eq,      0));

        // build the all-rows slack matrices:
        std::shared_ptr<sparse_matrix> jacobian_slacks_column0 = jacobian_slack_brow0_bcol0->vertical(jacobian_slack_brow1_bcol0);
        std::shared_ptr<sparse_matrix> jacobian_slacks_column1 = jacobian_slack_brow0_bcol1->vertical(jacobian_slack_brow1_bcol1);
        std::shared_ptr<sparse_matrix> jacobian_slacks_column2 = jacobian_slack_brow0_bcol2->vertical(jacobian_slack_brow1_bcol2);
        
        // horizontally-concatenate columns to make the full slack matrix:
        std::shared_ptr<sparse_matrix> jacobian_slacks_column0_column1 = jacobian_slacks_column0->horizontal(jacobian_slacks_column1);
        std::shared_ptr<sparse_matrix> jacobian_slacks_all             = jacobian_slacks_column0_column1->horizontal(jacobian_slacks_column2);
        
        // finally, build the full jacobian.
        jacobian_ = jacobian_eq_ieq->horizontal(jacobian_slacks_all);
        
        // TODO this is quite inefficient. at the *very* least, could keep jacobian_slacks_c0_c1_c2 around/static/something.
        // .... on the other hand, it wasn't really showing up in the profiling results before, so... who knows?
        // \note Bug killed here, Sept 24, 2013: wasn't making quite the correct jacobian slacks matrix.
        // @TODO decide if I like this... could make this way cleaner by not doing the left_{num rows} == right_{num rows} style checks.
        // But... then bugs become a lot lot lot harder to find.
        
        prior_constructed_jacobians_[iterate.get_serial()] = jacobian_;
        if (PRINT) std::cout << "jacobian_: " << jacobian_ << std::endl;
    }
    
    // LAGRANGIAN PART
        std::map<int, std::shared_ptr<matrix_base_class> >::iterator hessian_it = prior_constructed_hessians_.find(iterate.get_serial());
        if (hessian_it != prior_constructed_hessians_.end()) {hessian_ = hessian_it->second; std::cout << "HESS HIT\n";}
    else {
        // Here, we're constructing:
        //  /  H   0  \
        //  \  0  0*I /
        // Original H was num_nlp_variables_ x num_nlp_variables. Lower right corner has 2*n_eq + n_ieq slacks.
        // num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_
        std::shared_ptr<sparse_matrix>  hess_lower_left(new sparse_matrix(2*num_nlp_constraints_eq_ + num_nlp_constraints_ieq_, num_nlp_variables_, 0));    // null matrix for lower left
        std::shared_ptr<sparse_matrix> hess_upper_right(new sparse_matrix(num_nlp_variables_, 2*num_nlp_constraints_eq_ + num_nlp_constraints_ieq_, 0));  // null matrix for upper right.
        std::shared_ptr<sparse_matrix> hess_lower_right(sparse_matrix::diagonal_matrix(2*num_nlp_constraints_eq_ + num_nlp_constraints_ieq_, 2*num_nlp_constraints_eq_ + num_nlp_constraints_ieq_, 0.0));                      // 0.0*I_{2*num_qp_con} for lower right. (fixes bug in qpOASES...)
        // std::cout << "hess_right_bottom: " << hess_right_bottom << std::endl;
    
        std::shared_ptr<sparse_matrix> hess_left = nlp_hessian->vertical(hess_lower_left);
        // std::cout << "hess_left: " << hess_left << std::endl;
        std::shared_ptr<sparse_matrix> hess_right = hess_upper_right->vertical(hess_lower_right);
        // std::cout << "hess_right: " << hess_right << std::endl;
        hessian_ = hess_left->horizontal(hess_right);
        prior_constructed_hessians_[iterate.get_serial()] = hessian_;
    }
    // std::cout << "hessian_sparse_: " << hessian_sparse_ << std::endl;
    
}