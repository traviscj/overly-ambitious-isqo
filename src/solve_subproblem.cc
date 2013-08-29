#include <typeinfo>

#include "solve_subproblem.hh"

SolveQuadraticProgram::SolveQuadraticProgram(Nlp &nlp) : 
    FunctionWithNLPState(nlp), 
    example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual()), 
    first_(true)
    {}



iSQOStep SolveQuadraticProgram::operator()(iSQOSparseQuadraticSubproblem &sparse_subproblem) {
	int nWSR = 2000;
    operator_setup();

    qpOASES::SymSparseMat qpoases_hessian_(sparse_subproblem.num_qp_variables_, sparse_subproblem.num_qp_variables_, &sparse_subproblem.hessian_sparse_.row_indices_[0], &sparse_subproblem.hessian_sparse_.col_starts_[0], &sparse_subproblem.hessian_sparse_.vals_[0]);
    qpoases_hessian_.createDiagInfo();
    qpoases_hessian_.addToDiag(sparse_subproblem.hessian_shift_);
    qpOASES::SparseMatrix qpoases_jacobian_(sparse_subproblem.num_qp_constraints_, sparse_subproblem.num_qp_variables_, &sparse_subproblem.jacobian_sparse_.row_indices_[0], &sparse_subproblem.jacobian_sparse_.col_starts_[0], &sparse_subproblem.jacobian_sparse_.vals_[0]);
    
    qpOASES::returnValue ret;
	if (first_) {
		ret = example_.init( &qpoases_hessian_,
							 &sparse_subproblem.gradient_[0],
							 &qpoases_jacobian_,
							 &sparse_subproblem.lower_bound_[0],
							 &sparse_subproblem.upper_bound_[0],
							 &sparse_subproblem.jacobian_lower_bound_[0],
							 &sparse_subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	} else{
		ret = example_.hotstart(&qpoases_hessian_,
							 &sparse_subproblem.gradient_[0],
							 &qpoases_jacobian_,
							 &sparse_subproblem.lower_bound_[0],
							 &sparse_subproblem.upper_bound_[0],
							 &sparse_subproblem.jacobian_lower_bound_[0],
							 &sparse_subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	}
	return operator_finish(ret);
}

void SolveQuadraticProgram::operator_setup() {

	qpOASES::Options opt;

	opt.terminationTolerance = 1e-6;
	example_.setOptions(opt);
	example_.setPrintLevel(qpOASES::PL_NONE);
}
iSQOStep SolveQuadraticProgram::operator_finish(qpOASES::returnValue ret) {
	std::cerr.flush();
	std::cout.flush();
    
	size_t return_status = example_.getStatus();
	// getGlobalMessageHandler()->listAllMessages();
	if( ret != qpOASES::SUCCESSFUL_RETURN ){
        // sparse_subproblem.print();
        printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
		first_ = true;
		example_.reset();
	} else {
		if (first_) first_ = false;
	}
	iSQOStep step(nlp_->num_primal(), nlp_->num_dual_eq(), nlp_->num_dual_ieq());
	// std::cout << "address of step: " << &step << std::endl;

	// full_primal needs to hold primal variables and all slack variables
	std::vector<double> primal_and_slack_values(nlp_->num_primal() + 2*nlp_->num_dual());
	example_.getPrimalSolution( &primal_and_slack_values[0] );
	for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
		step.primal_values_[primal_index] = primal_and_slack_values[primal_index];

	// std::cout << "primal: " << step.primal_values_ << std::endl;
	// 	- every QP variable has a multiplier:
	//		-- iterate.num_primal_: primal variables in NLP
	// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
	//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
	//	- every constraint has a multiplier:
	//		-- iterate.num_dual_eq_
	//		-- iterate.num_dual_ieq_
	// so in total: 
	std::vector<double> yOpt(nlp_->num_primal() + 2*(nlp_->num_dual()) + nlp_->num_dual());
	example_.getDualSolution( &yOpt[0] );				

	step.status_ = ret;
	// to get eq con multipliers, read past NLP & slack variables:
	for (size_t dual_eq_index=0; dual_eq_index<nlp_->num_dual_eq(); ++dual_eq_index) {
		step.dual_eq_values_[dual_eq_index] = -yOpt[nlp_->num_primal() + 2*(nlp_->num_dual()) + dual_eq_index];
	}
	// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
	for (size_t dual_ieq_index=0; dual_ieq_index<nlp_->num_dual_ieq(); ++dual_ieq_index) {
		step.dual_ieq_values_[dual_ieq_index] = -yOpt[nlp_->num_primal() + 2*(nlp_->num_dual()) + nlp_->num_dual_eq() + dual_ieq_index];
	}

	bool PRINT=false;
	if (PRINT) {
		std::cout << step;
	}
	return step;
}

iSQOStep SolveQuadraticProgram::operator()(iSQOQuadraticSubproblem &subproblem) {
	int nWSR = 2000;

    qpOASES::SymDenseMat qpoases_hessian_(subproblem.num_qp_variables_, subproblem.num_qp_variables_, subproblem.num_qp_variables_, &subproblem.hessian_.data_[0]);
    qpOASES::DenseMatrix qpoases_jacobian_(subproblem.num_qp_constraints_, subproblem.num_qp_variables_, subproblem.num_qp_variables_, &subproblem.jacobian_.data_[0]);
	qpOASES::returnValue ret;

	if (first_) {
		ret = example_.init( &qpoases_hessian_,
							 &subproblem.gradient_[0],
							 &qpoases_jacobian_,
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	} else{
		ret = example_.hotstart(&qpoases_hessian_,
							 &subproblem.gradient_[0],
							 &qpoases_jacobian_,
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	}
	return operator_finish(ret);
}
