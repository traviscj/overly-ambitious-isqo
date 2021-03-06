// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)


#ifndef PENALTY_FUNCTION_HH_DSWMBAWM
#define PENALTY_FUNCTION_HH_DSWMBAWM

#include "nlp.hh"
#include "iterate.hh"
#include "constraint_violation.hh"

//! \brief a class returning a function to evaluate the penalty function.
//! 
//! The penalty function is:
//! \f[ \phi(x, \mu) := \mu f(x) + \| c_{\mathcal{E}}(x) \|_1 + \| \left[ c_{\mathcal{I}}(x) \right]^+ \|_1 \f]
class PenaltyFunction : public FunctionWithNLPState{
public:
	PenaltyFunction(iSQOControlPanel &control, Nlp &nlp);
	double operator()(const iSQOIterate &iterate) const;
protected:
	ConstraintViolationFunction constraint_violation_func_;
};


#endif /* end of include guard: PENALTY_FUNCTION_HH_DSWMBAWM */
