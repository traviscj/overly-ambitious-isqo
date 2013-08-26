
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
	return iterate.penalty_parameter_*f + constraint_violation_func_(iterate);
}