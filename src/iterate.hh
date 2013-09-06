
#ifndef __GUARD_ITERATE_HH
#define __GUARD_ITERATE_HH

#include <iostream>
#include <vector>

#include "step.hh"

//! \brief class for storing iSQO iterates (primal values and one set of dual values)
class iSQOIterate {
public:
	iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq);
	iSQOIterate(const iSQOIterate& s);
	double x_norm() const;
	std::ostream &print(std::ostream &) const;
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

//! \brief print an iterate
//!
//! here, we print out the underlying vectors associated with the iterate iter, to the output stream 'os'.
inline std::ostream& operator<< (std::ostream& os, const iSQOIterate& iter) {
    return iter.print(os);
}

#endif
