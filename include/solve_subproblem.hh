#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp);
	

    // iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
    
    qpOASES::SymmetricMatrix *get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> hessian);
    qpOASES::SymmetricMatrix *get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> hessian);
    qpOASES::SymmetricMatrix *get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> hessian);
    qpOASES::Matrix *get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> jacobian);
    qpOASES::Matrix *get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> jacobian);
    qpOASES::Matrix *get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> jacobian);
    
private:
    void operator_setup();
    iSQOStep operator_finish(const iSQOQuadraticSubproblem &subproblem, qpOASES::returnValue ret);
protected:
    // virtual void solve(iSQOQuadraticSubproblem *subproblem) = 0;
    
	qpOASES::SQProblem example_;
	bool first_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
