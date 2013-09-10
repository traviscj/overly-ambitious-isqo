#ifndef TEXT_OUTPUT_HH_45I5FVHZ
#define TEXT_OUTPUT_HH_45I5FVHZ

#include "step.hh"
#include "iterate.hh"
#include "subproblem.hh"
#include "nlp.hh"
#include "constraint_violation.hh"
#include "penalty_function.hh"
#include "linear_model_reduction.hh"
#include "residual_function.hh"
#include "nlp_state.hh"

//! \brief class responsible for doing text output as the algorithm progresses
class TextOutput : public FunctionWithNLPState {
public:
	TextOutput (Nlp &nlp);
    void nlp();
	void start();
	void pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const;
	void subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const;
	void subproblem_skip();
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, std::string step_type, double step_mix, double alpha);
protected:
	// Nlp *nlp_;
	ConstraintViolationFunction constraint_violation_func_;
	LinearReductionFunction linear_decrease_func_;
	PenaltyFunction pen_func_;
	ResidualFunction residual_func_;
	
	static const char output_desc_pre_[];
	static const char output_desc_subprob_[];
	static const char output_desc_post_[];
	static const char output_format_pre_[];
	static const char output_format_subprob_[];
	static const char output_format_post_[];
};


#endif /* end of include guard: TEXT_OUTPUT_HH_45I5FVHZ */
