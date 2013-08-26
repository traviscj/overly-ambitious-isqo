

#ifndef PENALTY_FUNCTION_HH_DSWMBAWM
#define PENALTY_FUNCTION_HH_DSWMBAWM

#include "nlp.hh"
#include "iterate.hh"
#include "constraint_violation.hh"

class PenaltyFunction : public FunctionWithNLPState{
public:
	PenaltyFunction(Nlp &nlp);
	double operator()(const iSQOIterate &iterate) const;
protected:
	ConstraintViolationFunction constraint_violation_func_;
};


#endif /* end of include guard: PENALTY_FUNCTION_HH_DSWMBAWM */
