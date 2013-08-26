
#ifndef __GUARD_CONSTRAINT_VIOLATION
#define __GUARD_CONSTRAINT_VIOLATION

#include <vector>

#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_state.hh"

class ConstraintViolationFunction : public FunctionWithNLPState {
public:
	ConstraintViolationFunction(Nlp &nlp);
	double operator()(const iSQOIterate &iterate) const;
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const;	
	double two_vectors(const std::vector<double> &eq_items, const std::vector<double> &ieq_items) const;
protected:
};

#endif
