// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

#define QPOASES_PROBLEM qpOASES::SQProblem
// #define QPOASES_PROBLEM qpOASES::SQProblemSchur

//! \brief class which returns a function for actually solving QPs with qpOASES
class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(iSQOControlPanel &control, Nlp &nlp);
	~SolveQuadraticProgram();

    // iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem, const iSQOIterate &iterate);

    void save_qp_state() {
        // backup_ = std::shared_ptr<qpOASES::SQProblemSchur>(new qpOASES::SQProblemSchur(*example_));
        backup_ = std::shared_ptr<QPOASES_PROBLEM>(new QPOASES_PROBLEM(*example_));
    }
    void restore_qp_state() {
        example_ = std::shared_ptr<QPOASES_PROBLEM>(new QPOASES_PROBLEM(*backup_));
    }
    
        
private:
    void operator_setup();
    iSQOStep operator_finish(const iSQOQuadraticSubproblem &subproblem, const iSQOIterate &iterate, int nWSR, qpOASES::returnValue ret);
    
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> hessian);
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> hessian);
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> hessian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> jacobian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> jacobian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> jacobian);
    
protected:
    // virtual void solve(iSQOQuadraticSubproblem *subproblem) = 0;
    
    std::shared_ptr<QPOASES_PROBLEM> example_;
    std::shared_ptr<QPOASES_PROBLEM> backup_;
    //  *example_;
    // QPOASES_PROBLEM *backup_;
    
    std::shared_ptr<qpOASES::Options> opt_;
  std::shared_ptr<qpOASES::SymSparseMat> qpoases_identity_for_hessian_;
  std::shared_ptr<qpOASES::Constraints> initial_constraints_;
  std::shared_ptr<qpOASES::Bounds> initial_bounds_;
  //std::vector<long int> qpoases_hessian_diag_info_;
  int *qpoases_hessian_diag_info_;
	bool first_;
    bool last_successful_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
