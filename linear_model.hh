#ifndef LINEAR_MODEL_HH_42HS2XKK
#define LINEAR_MODEL_HH_42HS2XKK

#include "iterate.hh"
#include "step.hh"
#include "nlp.hh"
#include "nlp_state.hh"
#include "constraint_violation.hh"

class LinearModelFunction : public FunctionWithNLPState {
public:
	LinearModelFunction(Nlp &nlp);
	double operator()(const iSQOIterate &iterate, const iSQOStep &step);
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

#endif /* end of include guard: LINEAR_MODEL_HH_42HS2XKK */
