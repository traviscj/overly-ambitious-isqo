#ifndef RESIDUAL_FUNCTION_HH_QBT6Z7GO
#define RESIDUAL_FUNCTION_HH_QBT6Z7GO

#include <vector>

#include "iterate.hh"
#include "nlp.hh"
#include "subproblem.hh"
#include "step.hh"

class ResidualFunction : public FunctionWithNLPState {
public:
	ResidualFunction(Nlp &nlp);
	double operator()(const iSQOIterate &iterate) const;
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const;
protected:
	double resid_helper(const iSQOIterate &iterate, std::vector<double> stationarity, std::vector<double> constraint_eq_values, std::vector<double> constraint_ieq_values, std::vector<double> constraint_eq_dual_values, std::vector<double> constraint_ieq_dual_values) const;
    
};

#endif /* end of include guard: RESIDUAL_FUNCTION_HH_QBT6Z7GO */
