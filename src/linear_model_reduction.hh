#ifndef LINEAR_MODEL_REDUCTION_HH_CSYGSX4Y
#define LINEAR_MODEL_REDUCTION_HH_CSYGSX4Y

#include "iterate.hh"
#include "step.hh"
#include "nlp.hh"
#include "nlp_state.hh"
#include "constraint_violation.hh"

//! \brief class returning a function to evaluate linear model reduction
class LinearReductionFunction : public FunctionWithNLPState {
public:
	LinearReductionFunction(Nlp &nlp);
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const;
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

#endif /* end of include guard: LINEAR_MODEL_REDUCTION_HH_CSYGSX4Y */
