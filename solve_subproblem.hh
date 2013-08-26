#ifndef SOLVE_SUBPROBLEM_HH_9BLKFD95
#define SOLVE_SUBPROBLEM_HH_9BLKFD95

#include <qpOASES.hpp>

#include "subproblem.hh"
#include "nlp_state.hh"

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp);
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem);
private:
protected:
	qpOASES::SQProblem example_;
	bool first_;
};

#endif /* end of include guard: SOLVE_SUBPROBLEM_HH_9BLKFD95 */
