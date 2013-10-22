// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#include "penalty_function.hh"

PenaltyFunction::PenaltyFunction(Nlp &nlp) : 
		FunctionWithNLPState(nlp), 
		constraint_violation_func_(nlp)
		{
	// std::cout << "-- Initializing penalty function." << std::endl;
}
double PenaltyFunction::operator()(const iSQOIterate &iterate) const {
	// std::cout << "--- Calling the PenaltyFunction Functor..." << std::endl;
	double f = nlp_->objective(iterate);
	// std::cout << "--- objective: " << f << std::endl;
	return iterate.get_penalty_parameter()*f + constraint_violation_func_(iterate);
}