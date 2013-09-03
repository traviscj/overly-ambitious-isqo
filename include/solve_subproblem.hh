#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp);
	

	iSQOStep operator()(iSQOSparseQuadraticSubproblem &subproblem);
	iSQOStep operator()(iSQOQuadraticSubproblem &subproblem);
private:
    void operator_setup();
    iSQOStep operator_finish(const iSQOQuadraticSubproblem &subproblem, qpOASES::returnValue ret);
protected:
    // virtual void solve(iSQOQuadraticSubproblem *subproblem) = 0;
    
	qpOASES::SQProblem example_;
	bool first_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
