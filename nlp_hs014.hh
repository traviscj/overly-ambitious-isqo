
#ifndef __GUARD_NLP_HS014_HH
#define __GUARD_NLP_HS014_HH
#include <vector>
#include "iterate.hh"

class Hs014 : public Nlp {
public:
	Hs014();	
	iSQOIterate initial();
	double objective(const iSQOIterate &iterate);
	std::vector<double> constraints_equality(const iSQOIterate &iterate);
	std::vector<double> constraints_inequality(const iSQOIterate &iterate);
	
	std::vector<double> objective_gradient(const iSQOIterate &iterate);
	matrix constraints_equality_jacobian(const iSQOIterate &iterate);
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate);
	
	matrix lagrangian_hessian(const iSQOIterate &iterate);
protected:
};
#endif
