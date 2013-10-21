
#include "penalty_function.hh"

PenaltyFunction::PenaltyFunction(iSQOControlPanel &control, Nlp &nlp) : 
		FunctionWithNLPState(control, nlp), 
		constraint_violation_func_(control, nlp)
		{
	// std::cout << "-- Initializing penalty function." << std::endl;
}
double PenaltyFunction::operator()(const iSQOIterate &iterate) const {
	// std::cout << "--- Calling the PenaltyFunction Functor..." << std::endl;
	double f = nlp_->objective(iterate);
	// std::cout << "--- objective: " << f << std::endl;
	return iterate.get_penalty_parameter()*f + constraint_violation_func_(iterate);
}