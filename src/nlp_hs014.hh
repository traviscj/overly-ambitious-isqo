// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#ifndef __GUARD_NLP_HS014_HH
#define __GUARD_NLP_HS014_HH
#include <vector>
#include "iterate.hh"

//! \brief a specific Nlp class for the problem Hs014, for demonstration & testing.
class Hs014 : public Nlp {
public:
	Hs014();	
	iSQOIterate initial();
	double objective(const iSQOIterate &iterate);
	std::vector<double> constraints_equality(const iSQOIterate &iterate);
	std::vector<double> constraints_inequality(const iSQOIterate &iterate);
	
	std::vector<double> objective_gradient(const iSQOIterate &iterate);
	std::shared_ptr<matrix_base_class> constraints_equality_jacobian(const iSQOIterate &iterate);
	std::shared_ptr<matrix_base_class> constraints_inequality_jacobian(const iSQOIterate &iterate);
	std::shared_ptr<matrix_base_class> lagrangian_hessian(const iSQOIterate &iterate);
protected:
};
#endif
