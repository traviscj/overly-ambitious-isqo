// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)
#ifndef LINE_SEARCH_HH_KXP2YSCA
#define LINE_SEARCH_HH_KXP2YSCA

#include "step.hh"
#include "iterate.hh"
#include "nlp.hh"
#include "nlp_state.hh"
#include "penalty_function.hh"
#include "linear_model_reduction.hh"

//! \brief class returning a function to perform line search
class LineSearchFunction : public FunctionWithNLPState {
public:
	LineSearchFunction(iSQOControlPanel &control, Nlp &nlp);
	double operator()(const iSQOIterate &iterate, const iSQOStep &step);
private:
protected:
	PenaltyFunction penalty_func_;
	LinearReductionFunction linear_decrease_func_;
};

#endif /* end of include guard: LINE_SEARCH_HH_KXP2YSCA */
