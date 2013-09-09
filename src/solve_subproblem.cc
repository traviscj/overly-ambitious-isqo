#include <memory>
#include <typeinfo>

#include "solve_subproblem.hh"

SolveQuadraticProgram::SolveQuadraticProgram(Nlp &nlp) : 
    FunctionWithNLPState(nlp), 
    example_(0),
    first_(true)
    {
        // TODO: switch this based on sparsity.
        // TODO: smart pointer or deconstruct this.
        example_ = std::shared_ptr<qpOASES::SQProblem>(new qpOASES::SQProblemSchur(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual(), qpOASES::HST_UNKNOWN));
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
    
    return std::shared_ptr<qpOASES::SymmetricMatrix>(new qpOASES::SymDenseMat(subproblem.num_qp_variables_, subproblem.num_qp_variables_, subproblem.num_qp_variables_, &(*hessian->get_data())[0]));
}

std::shared_ptr<qpOASES::SymmetricMatrix> SolveQuadraticProgram::get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> hessian) {
    // TODO fix this casting, but later...
    std::shared_ptr<qpOASES::SymSparseMat> qpoases_hessian(new qpOASES::SymSparseMat(   
                                                                subproblem.num_qp_variables_, 
                                                                subproblem.num_qp_variables_, 
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
    return std::shared_ptr<qpOASES::Matrix>(new qpOASES::DenseMatrix(subproblem.num_qp_constraints_, subproblem.num_qp_variables_, subproblem.num_qp_variables_, &(*jacobian->get_data())[0]));
}

std::shared_ptr<qpOASES::Matrix> SolveQuadraticProgram::get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> jacobian) {
    // need to cheat here! don't mess with my matrix, qpOASES!
    std::shared_ptr<sparse_matrix> nonconst_jacobian(std::const_pointer_cast<sparse_matrix>(jacobian));
    std::shared_ptr<qpOASES::Matrix> qpoases_jacobian(new qpOASES::SparseMatrix(
                                                                subproblem.num_qp_constraints_, 
                                                                subproblem.num_qp_variables_, 
                                                                (int*)&nonconst_jacobian->get_row_indices()[0],
                                                                (int*)&nonconst_jacobian->get_col_starts()[0], 
                                                                (double*)&nonconst_jacobian->get_vals()[0]
                                                                    ));
    return qpoases_jacobian;
}

// DENSE quadratic subproblem solver:
iSQOStep SolveQuadraticProgram::operator()(iSQOQuadraticSubproblem &subproblem) {
	int nWSR = 50000;
    operator_setup();
    
    // std::cout << std::endl << "immediately before solve: " << subproblem << std::endl;
    std::shared_ptr<qpOASES::SymmetricMatrix> qpoases_hessian_ = this->get_qpoases_hessian(subproblem, subproblem.hessian_);
    std::shared_ptr<qpOASES::Matrix> qpoases_jacobian_ = this->get_qpoases_jacobian(subproblem, subproblem.jacobian_);
	qpOASES::returnValue ret;
// subproblem.print();
	if (first_) {
		ret = example_->init( qpoases_hessian_.get(),
							 &subproblem.gradient_[0],
							 qpoases_jacobian_.get(),
							 &subproblem.lower_bound_[0],
							 &subproblem.upper_bound_[0],
							 &subproblem.jacobian_lower_bound_[0],
							 &subproblem.jacobian_upper_bound_[0],
							 nWSR,
							 0);
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
	return operator_finish(subproblem, ret);
}

// SPARSE quadratic subproblem solver:
// iSQOStep SolveQuadraticProgram::operator()(iSQOSparseQuadraticSubproblem &sparse_subproblem) {
//     int nWSR = 2000;
//     operator_setup();
// 
//     qpOASES::SymSparseMat qpoases_hessian_();
//     
//     qpOASES::SparseMatrix qpoases_jacobian_(sparse_subproblem.num_qp_constraints_, sparse_subproblem.num_qp_variables_, &sparse_subproblem.jacobian_sparse_.row_indices_[0], &sparse_subproblem.jacobian_sparse_.col_starts_[0], &sparse_subproblem.jacobian_sparse_.vals_[0]);
//     
//     qpOASES::returnValue ret;
//     if (first_) {
//         ret = example_->init( &qpoases_hessian_,
//                              &sparse_subproblem.gradient_[0],
//                              &qpoases_jacobian_,
//                              &sparse_subproblem.lower_bound_[0],
//                              &sparse_subproblem.upper_bound_[0],
//                              &sparse_subproblem.jacobian_lower_bound_[0],
//                              &sparse_subproblem.jacobian_upper_bound_[0],
//                              nWSR,
//                              0);
//     } else{
//         ret = example_->hotstart(&qpoases_hessian_,
//                              &sparse_subproblem.gradient_[0],
//                              &qpoases_jacobian_,
//                              &sparse_subproblem.lower_bound_[0],
//                              &sparse_subproblem.upper_bound_[0],
//                              &sparse_subproblem.jacobian_lower_bound_[0],
//                              &sparse_subproblem.jacobian_upper_bound_[0],
//                              nWSR,
//                              0);
//     }
//     return operator_finish(sparse_subproblem, ret);
// }

void SolveQuadraticProgram::operator_setup() {

    // qpOASES::Options *opt = new qpOASES::Options;
    // 
    //     // opt.terminationTolerance = 1e-6;
    //     opt->setToReliable();
    //     
    //     opt->enableEqualities = qpOASES::BooleanType(true);
    //     opt->enableFullLITests = qpOASES::BooleanType(true);
    //     // opt->epsLITests ; // maybe play with this later. (maybe 1e-8)
    //     opt->enableRegularisation = qpOASES::BooleanType(true);
    //     opt->enableNZCTests = qpOASES::BooleanType(true);
    //     // opt->initialStatusBounds = qpOASES::ST_LOWER;
    //     // opt->numRegularisationSteps = 200;
    //     // std::cout << std::endl << "max reg steps: " <<  opt.numRegularisationSteps << std::endl;
    //     example_->setOptions(*opt);
	example_->setPrintLevel(qpOASES::PL_NONE);
}
iSQOStep SolveQuadraticProgram::operator_finish(const iSQOQuadraticSubproblem &subproblem, qpOASES::returnValue ret) {
	std::cerr.flush();
	std::cout.flush();
    
	size_t return_status = example_->getStatus();
	// getGlobalMessageHandler()->listAllMessages();
	if( ret != qpOASES::SUCCESSFUL_RETURN ){
        // 
        if (ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS || 
            ret == qpOASES::RET_QP_UNBOUNDED || 
            ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS || 
            ret == qpOASES::RET_STEPDIRECTION_FAILED_CHOLESKY ||      // only in dense.
            ret == qpOASES::RET_INIT_FAILED_HOTSTART
            ) {
            // don't need to do anything for unbounded QPs...
             
        } else {
            // subproblem.print();
            // 
            std::cout << "failed subproblem: " << subproblem << std::endl;
            printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
        
        }
        first_ = true;
        example_->reset();
	} else {
		if (first_) first_ = false;
	}
	iSQOStep step(nlp_->num_primal(), nlp_->num_dual_eq(), nlp_->num_dual_ieq());
	// std::cout << "address of step: " << &step << std::endl;

	// full_primal needs to hold primal variables and all slack variables
	std::vector<double> primal_and_slack_values(nlp_->num_primal() + 2*nlp_->num_dual());
	example_->getPrimalSolution( &primal_and_slack_values[0] );
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
	// so in total:          VARIABLE MULTIPLIERS   SLACK MULTIPLIERS    CONSTRAINT MULTIPLIERS
	std::vector<double> yOpt(nlp_->num_primal() + 2*(nlp_->num_dual()) + nlp_->num_dual());
	example_->getDualSolution( &yOpt[0] );				

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
    
    // assert(false);
	return step;
}

