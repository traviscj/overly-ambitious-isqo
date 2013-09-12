#include <memory>
#include <typeinfo>

#include "solve_subproblem.hh"

SolveQuadraticProgram::SolveQuadraticProgram(Nlp &nlp) : 
    FunctionWithNLPState(nlp), 
    example_(NULL),
    backup_(NULL),
    opt_(NULL),
    first_(true),
    last_successful_(false)
    {
        // TODO: switch this based on sparsity.
        // TODO: smart pointer or deconstruct this.
        int num_qp_con = (int)(nlp_->num_dual());
        int num_qp_var = (int)(nlp_->num_primal() + 2*nlp_->num_dual());
        
        example_ = std::shared_ptr<qpOASES::SQProblemSchur>(new qpOASES::SQProblemSchur(num_qp_var, num_qp_con, qpOASES::HST_UNKNOWN));
        opt_ = std::shared_ptr<qpOASES::Options>(new qpOASES::Options);
        // 
        //     // opt.terminationTolerance = 1e-6;
        // opt->setToReliable();
        opt_->enableWorkingSetRepair = qpOASES::BooleanType(true);
        opt_->enableEqualities = qpOASES::BooleanType(true);
        // opt_->enableFarBounds = qpOASES::BooleanType(false);
        opt_->enableFullLITests = qpOASES::BooleanType(true);
        opt_->epsLITests = 1e-12; // maybe play with this later. (maybe 1e-8)
        // usually too low, rather than too high, but might help with the 'infeasible QP' troubles.
        //     opt->enableRegularisation = qpOASES::BooleanType(true);
        opt_->enableNZCTests = qpOASES::BooleanType(true);
        // opt_->initialStatusBounds = qpOASES::ST_INACTIVE;
        // opt_->initialStatusBounds = qpOASES::ST_LOWER;
        opt_->printLevel = qpOASES::PL_NONE;
        // opt_->printLevel = qpOASES::PL_MEDIUM;
    
        //     // opt->numRegularisationSteps = 200;
        //     // std::cout << std::endl << "max reg steps: " <<  opt.numRegularisationSteps << std::endl;
        example_->setOptions(*opt_);
        // example_->setPrintLevel(qpOASES::PL_NONE);
            // example_->setPrintLevel(qpOASES::PL_MEDIUM);
        
    }

SolveQuadraticProgram::~SolveQuadraticProgram() {
    
}

std::shared_ptr<qpOASES::SymmetricMatrix> SolveQuadraticProgram::get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> hessian) {
    std::shared_ptr<dense_matrix> attempt_dense = std::dynamic_pointer_cast< dense_matrix >(hessian);
    std::shared_ptr<sparse_matrix> attempt_sparse = std::dynamic_pointer_cast< sparse_matrix >(hessian);
    
    if (attempt_dense != NULL) {
        return get_qpoases_hessian(subproblem, attempt_dense);
    } else if (attempt_sparse != NULL) {
        return get_qpoases_hessian(subproblem, attempt_sparse);
    } else {
        throw("wtf mate");
    }
}
std::shared_ptr<qpOASES::SymmetricMatrix> SolveQuadraticProgram::get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> hessian) {
    // std::vector<double> *data_vec_double = ;
    // double *data_double_star = ;
    // int num_qp_con = (int)subproblem.num_qp_constraints_;
    int num_qp_var = (int)subproblem.num_qp_variables_;
    
    return std::shared_ptr<qpOASES::SymmetricMatrix>(new qpOASES::SymDenseMat(num_qp_var, num_qp_var, num_qp_var, &(*hessian->get_data())[0]));
}

std::shared_ptr<qpOASES::SymmetricMatrix> SolveQuadraticProgram::get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> hessian) {
    // TODO fix this casting, but later...
    // int num_qp_con = (int)subproblem.num_qp_constraints_;
    // std::cout << "SolveQuadraticProgram::get_qpoases_hessian: sizeof subproblem.num_qp_constraints: " << sizeof(subproblem.num_qp_constraints_);
    // std::cout << "; sizeof num_qp_con: " << sizeof(num_qp_con) << std::endl;
    
    int num_qp_var = (int)subproblem.num_qp_variables_;
    
    std::shared_ptr<qpOASES::SymSparseMat> qpoases_hessian(new qpOASES::SymSparseMat(   
                                                                num_qp_var, 
                                                                num_qp_var, 
                                                                (int*)(&hessian->get_row_indices()[0]),
                                                                (int*)(&hessian->get_col_starts()[0]), 
                                                                (double*)(&hessian->get_vals()[0])
                                                                    ));
    qpoases_hessian->createDiagInfo();
    return qpoases_hessian;
    // return (new qpOASES::SymSparseMat(subproblem.num_qp_variables_, subproblem.num_qp_variables_, subproblem.num_qp_variables_, &hessian->data_[0]));
}
std::shared_ptr<qpOASES::Matrix> SolveQuadraticProgram::get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> jacobian) {
    // std::shared_ptr<dense_matrix> attempt_dense = std::dynamic_pointer_cast< dense_matrix  >(nlp_eq_jacobian);
    // TODO fix this, before anyone notices!
    std::shared_ptr<dense_matrix> attempt_dense = std::dynamic_pointer_cast< dense_matrix >(jacobian);
    std::shared_ptr<sparse_matrix> attempt_sparse = std::dynamic_pointer_cast< sparse_matrix >(jacobian);
    
    if (attempt_dense != NULL) {
        return get_qpoases_jacobian(subproblem, attempt_dense);
    } else if (attempt_sparse != NULL) {
        return get_qpoases_jacobian(subproblem, attempt_sparse);
    } else {
        throw("wtf mate");
    }
}
std::shared_ptr<qpOASES::Matrix> SolveQuadraticProgram::get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> jacobian) {
    int num_qp_con = subproblem.num_qp_constraints_;
    int num_qp_var = subproblem.num_qp_variables_;
    return std::shared_ptr<qpOASES::Matrix>(new qpOASES::DenseMatrix(num_qp_con, num_qp_var, subproblem.num_qp_variables_, &(*jacobian->get_data())[0]));
}

std::shared_ptr<qpOASES::Matrix> SolveQuadraticProgram::get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> jacobian) {
    // need to cheat here! don't mess with my matrix, qpOASES!
    std::shared_ptr<sparse_matrix> nonconst_jacobian(std::const_pointer_cast<sparse_matrix>(jacobian));
    int num_qp_con = (int)subproblem.num_qp_constraints_;
        
    
    int num_qp_var = (int)subproblem.num_qp_variables_;
    std::shared_ptr<qpOASES::Matrix> qpoases_jacobian(new qpOASES::SparseMatrix(
                                                                num_qp_con, 
                                                                num_qp_var, 
                                                                (int*)&nonconst_jacobian->get_row_indices()[0],
                                                                (int*)&nonconst_jacobian->get_col_starts()[0], 
                                                                (double*)&nonconst_jacobian->get_vals()[0]
                                                                    ));
    return qpoases_jacobian;
}

// DENSE quadratic subproblem solver:
iSQOStep SolveQuadraticProgram::operator()(iSQOQuadraticSubproblem &subproblem) {
	int nWSR = 50000;
    if (!first_ && !last_successful_) {
        restore_qp_state();
    }
    operator_setup();
    
    // std::cout << std::endl << "immediately before solve: " << subproblem << std::endl;
    std::shared_ptr<qpOASES::SymmetricMatrix> qpoases_hessian_ = this->get_qpoases_hessian(subproblem, subproblem.hessian_);
    std::shared_ptr<qpOASES::Matrix> qpoases_jacobian_ = this->get_qpoases_jacobian(subproblem, subproblem.jacobian_);
	qpOASES::returnValue ret;
    example_->setOptions(*opt_);
    double cputime=0;
// subproblem.print();
	if (first_) {
        
        double *initial_guess_primal = new double[nlp_->num_primal() + 2*nlp_->num_dual()];
        double *initial_guess_dual = new double[nlp_->num_dual()];
        qpOASES::Bounds *initial_bounds = new qpOASES::Bounds((int)(nlp_->num_primal() + 2*nlp_->num_dual()));
        initial_bounds->setupAllFree();
        for (int variable_index=(int)(nlp_->num_primal()); variable_index < (int)(nlp_->num_primal() + 2*nlp_->num_dual()); ++variable_index)
            initial_bounds->setStatus(variable_index, qpOASES::ST_LOWER);
        qpOASES::Constraints *initial_constraints = new qpOASES::Constraints((int)nlp_->num_dual());
        initial_constraints->setupAllLower();
        // for (int constraint_index=0; constraint_index < (int)(nlp_->num_dual()); ++constraint_index)
            // initial_constraints->setStatus(constraint_index, qpOASES::ST_LOWER);
        
		ret = example_->init( qpoases_hessian_.get(),
							 &subproblem.gradient_[0],
							 qpoases_jacobian_.get(),
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
                             &cputime,
                             0,
                             0,
                             initial_bounds,
                             initial_constraints);
	} else{
        
		ret = example_->hotstart(qpoases_hessian_.get(),
							 &subproblem.gradient_[0],
							 qpoases_jacobian_.get(),
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
	}
	return operator_finish(subproblem, nWSR, ret);
}

void SolveQuadraticProgram::operator_setup() {
}
iSQOStep SolveQuadraticProgram::operator_finish(const iSQOQuadraticSubproblem &subproblem, int nWSR, qpOASES::returnValue ret) {
	std::cerr.flush();
	std::cout.flush();
    
    // int return_status = example_->getStatus();
	// getGlobalMessageHandler()->listAllMessages();
	if( ret != qpOASES::SUCCESSFUL_RETURN ){
        // 
        if (( ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS || 
            ret == qpOASES::RET_QP_UNBOUNDED || 
            ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS || 
            ret == qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY ||      // only in dense.
            ret == qpOASES::RET_INIT_FAILED_HOTSTART
            )) {
            // don't need to do anything for unbounded QPs...
             
        } else {
            std::cout << "failed subproblem: " << subproblem << std::endl;
            printf( "# Working Set: [%d]\n", nWSR );
            printf( "Decoded Error Message: [%s]\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
        
        }
        if (first_) example_->reset();
        // first_ = false;
        last_successful_ = false;
	} else {
        first_ = false;
		last_successful_ = true;
	}
    
	iSQOStep step(nlp_->num_primal(), nlp_->num_dual_eq(), nlp_->num_dual_ieq(), ret, nWSR);
	// std::cout << "address of step: " << &step << std::endl;

	// full_primal needs to hold primal variables and all slack variables
	std::vector<double> primal_and_slack_values(nlp_->num_primal() + 2*nlp_->num_dual());
	example_->getPrimalSolution( &primal_and_slack_values[0] );
	for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
		step.set_primal_value(primal_index, primal_and_slack_values[primal_index]);

	// std::cout << "primal: " << step.primal_values_ << std::endl;
	// 	- every QP variable has a multiplier:
	//		-- iterate.num_primal_: primal variables in NLP
	// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
	//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
	//	- every constraint has a multiplier:
	//		-- iterate.num_dual_eq_
	//		-- iterate.num_dual_ieq_
	// so in total:          VARIABLE MULTIPLIERS   SLACK MULTIPLIERS    CONSTRAINT MULTIPLIERS
	std::vector<double> yOpt(nlp_->num_primal() + 2*(nlp_->num_dual()) + nlp_->num_dual());
	example_->getDualSolution( &yOpt[0] );				

	// to get eq con multipliers, read past NLP & slack variables:
	for (size_t dual_eq_index=0; dual_eq_index<nlp_->num_dual_eq(); ++dual_eq_index) {
		step.set_dual_eq_value(dual_eq_index, -yOpt[nlp_->num_primal() + 2*(nlp_->num_dual()) + dual_eq_index]);
	}
	// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
	for (size_t dual_ieq_index=0; dual_ieq_index<nlp_->num_dual_ieq(); ++dual_ieq_index) {
		step.set_dual_ieq_value(dual_ieq_index, -yOpt[nlp_->num_primal() + 2*(nlp_->num_dual()) + nlp_->num_dual_eq() + dual_ieq_index]);
	}

	bool PRINT=false;
	if (PRINT) {
		std::cout << step;
	}
    
    // assert(false);
	return step;
}
