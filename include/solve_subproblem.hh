#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp);
	
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
private:
protected:
    // virtual void solve(iSQOQuadraticSubproblem *subproblem) = 0;
    
	qpOASES::SQProblem example_;
	bool first_;
};

class SolveDenseQuadraticProgram : public SolveQuadraticProgram {
public:
	SolveDenseQuadraticProgram(Nlp &nlp);
    // iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
private:
protected:
    void solve(iSQOQuadraticSubproblem &subproblem);
    // qpOASES::SymDenseMat *qpoases_hessian_;
    // qpOASES::DenseMatrix *qpoases_jacobian_;
};

class SolveSparseQuadraticProgram : public SolveQuadraticProgram {
public:
	SolveSparseQuadraticProgram(Nlp &nlp);
    // iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
private:
protected:
    void solve(iSQOSparseQuadraticSubproblem &subproblem);
    // qpOASES::SymSparseMat *qpoases_hessian_;
    // qpOASES::SparseMatrix *qpoases_jacobian_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
