
#ifndef __GUARD_ITERATE_HH
#define __GUARD_ITERATE_HH

#include <iostream>
#include <vector>

#include "step.hh"

class iSQOIterate {
public:
	iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq);
	iSQOIterate(const iSQOIterate& s);
	double x_norm() const;
	std::ostream &print(std::ostream &);
	void update(const iSQOIterate &iterate, double alpha, const iSQOStep& step);
	void update_dual(const iSQOStep& step);
    
// private:
// protected:
	size_t num_primal_;
	size_t num_dual_eq_;
	size_t num_dual_ieq_;
	
	std::vector<double> primal_values_;
	std::vector<double> dual_eq_values_;
	std::vector<double> dual_ieq_values_;
	double penalty_parameter_;
};

#endif
