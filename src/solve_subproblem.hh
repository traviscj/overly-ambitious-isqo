#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

//! \brief class which returns a function for actually solving QPs with qpOASES
class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp);
	~SolveQuadraticProgram();

    // iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
        
private:
    void operator_setup();
    iSQOStep operator_finish(const iSQOQuadraticSubproblem &subproblem, qpOASES::returnValue ret);
    
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> hessian);
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> hessian);
    std::shared_ptr<qpOASES::SymmetricMatrix> get_qpoases_hessian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> hessian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<matrix_base_class> jacobian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<dense_matrix> jacobian);
    std::shared_ptr<qpOASES::Matrix> get_qpoases_jacobian(iSQOQuadraticSubproblem &subproblem, std::shared_ptr<sparse_matrix> jacobian);
    
protected:
    // virtual void solve(iSQOQuadraticSubproblem *subproblem) = 0;
    
	std::shared_ptr<qpOASES::SQProblem> example_;
	bool first_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
