
#include "solve_subproblem.hh"

SolveQuadraticProgram::SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual()), first_(true)  {}

iSQOStep SolveQuadraticProgram::operator()(const iSQOQuadraticSubproblem &subproblem) {
	// qpOASES::SQProblem example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual());
	/* Setting up QProblem object. */
	int nWSR = 2000;
	qpOASES::returnValue ret;
	// std::cout << std::endl;
	// std::cout << "hessian: " << subproblem.hessian_ <<std::endl;
	// std::cout << "jacobian: " << subproblem.jacobian_ <<std::endl;
	// std::cout << "gradient: " << subproblem.gradient_ << std::endl;
	// std::cout << "lb: : " << subproblem.lower_bound_ << std::endl;
	// std::cout << "ub: : " << subproblem.upper_bound_ << std::endl;
	// std::cout << "lbA: " << subproblem.jacobian_lower_bound_ << std::endl;
	// std::cout << "ubA: " << subproblem.jacobian_upper_bound_ << std::endl;
	
	qpOASES::Options opt;
	// opt.setToReliable();
	// opt.enableRegularisation = qpOASES::BooleanType(true);

	opt.terminationTolerance = 1e-6;
	example_.setOptions(opt);
	example_.setPrintLevel(qpOASES::PL_NONE);
	// std::cout.flush();
	// std::cerr.flush();
	if (first_) {
		// std::cout << "subproblem.hessian_.data_[0]: " << subproblem.hessian_ << std::endl;
		// std::cout << "subproblem: " << subproblem.gradient_ << std::endl;
		// std::cout << "subproblem: " << subproblem.jacobian_ << std::endl;
		// std::cout << "subproblem: " << subproblem.lower_bound_ << std::endl;
		// std::cout << "subproblem: " << subproblem.upper_bound_ << std::endl;
		// std::cout << "subproblem: " << subproblem.jacobian_lower_bound_ << std::endl;
		// std::cout << "subproblem: " << subproblem.jacobian_upper_bound_ << std::endl;
		// subproblem.print();
		ret = example_.init( &subproblem.hessian_.data_[0],
							 &subproblem.gradient_[0],
							 &subproblem.jacobian_.data_[0],
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	} else{
		
		// subproblem.print();
		ret = example_.hotstart(&subproblem.hessian_.data_[0],
							 &subproblem.gradient_[0],
							 &subproblem.jacobian_.data_[0],
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	}
	std::cerr.flush();
	std::cout.flush();

	size_t return_status = example_.getStatus();
	// getGlobalMessageHandler()->listAllMessages();
	if( ret != qpOASES::SUCCESSFUL_RETURN ){
		// subproblem.print();
        // printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
		first_ = true;
		example_.reset();
	} else {
		if (first_) first_ = false;
	}
	iSQOStep step(subproblem.num_nlp_variables_,subproblem.num_nlp_constraints_eq_,subproblem.num_nlp_constraints_ieq_);
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
	std::vector<double> yOpt(subproblem.num_qp_variables_ + subproblem.num_qp_constraints_);
	example_.getDualSolution( &yOpt[0] );				
	
	step.status_ = ret;
	// to get eq con multipliers, read past NLP & slack variables:
	for (size_t dual_eq_index=0; dual_eq_index<subproblem.num_nlp_constraints_eq_; ++dual_eq_index) {
		step.dual_eq_values_[dual_eq_index] = -yOpt[subproblem.num_qp_variables_+dual_eq_index];
	}
	// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
	for (size_t dual_ieq_index=0; dual_ieq_index<subproblem.num_nlp_constraints_ieq_; ++dual_ieq_index) {
		step.dual_ieq_values_[dual_ieq_index] = -yOpt[subproblem.num_qp_variables_+subproblem.num_nlp_constraints_eq_+dual_ieq_index];
	}
	// std::cout << "dual eq: " << step.dual_eq_values_ << std::endl;
	// std::cout << "dual ieq: " << step.dual_ieq_values_ << std::endl;
	
	bool PRINT=false;
	if (PRINT) {
		std::cout << step;
	}
	return step;
}
