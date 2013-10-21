
#ifndef __GUARD_ITERATE_HH
#define __GUARD_ITERATE_HH

#include <iostream>
#include <vector>

#include "step.hh"

static size_t global_iSQOIterate_serial = 0;

//! \brief class for storing iSQO iterates (primal values and one set of dual values)
class iSQOIterate {
public:
	iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq, double penalty_paramter);
//     iSQOIterate(const std::vector<double> &primals, const std::vector<double> &dual_eq, const std::vector<double> &dual_ieq, double penalty_parameter) : num_primal_(primals.size()),
// num_dual_eq_(dual_eq.size()),
// num_dual_ieq_(dual_ieq.size()),
// primal_values_(num_primal_),
// dual_eq_values_(num_dual_eq_),
// dual_ieq_values_(num_dual_ieq_),
// penalty_parameter_(penalty_parameter)
//  {
//         // primal_values_.assign(&primals[0], &primals[0] + num_primal_);
//         // dual_eq_values_.assign(&dual_eq[0], &dual_eq[0] + num_dual_eq_);
//         // dual_ieq_values_.assign(&dual_ieq[0], &dual_ieq[0] + num_dual_ieq_);
//          // assign_primal(&primals[0], &primals[0] + num_primal_);
//          // assign_dual_eq(&dual_eq[0], &dual_eq[0] + num_dual_eq_);
//          // assign_dual_ieq(&dual_ieq[0], &dual_ieq[0] + num_dual_ieq_);
//      assign_primal(primals);
//      assign_dual_eq(dual_eq);
//      assign_dual_ieq(dual_ieq);
//     }
	iSQOIterate(const iSQOIterate& s);
	double x_norm() const;
	std::ostream &print(std::ostream &) const;
	void update(const iSQOIterate &iterate, double alpha, const iSQOStep& step);
	void update_dual(const iSQOStep& step);
    inline double get_penalty_parameter() const { return penalty_parameter_; }
    inline const std::vector<double> *get_primal_values()   const { return &primal_values_;   }
    inline const std::vector<double> *get_dual_eq_values()  const { return &dual_eq_values_;  }
    inline const std::vector<double> *get_dual_ieq_values() const { return &dual_ieq_values_; }
    
    void set_penalty_parameter(double new_penalty_parameter) { penalty_parameter_ = new_penalty_parameter; }
    
    void assign_primal(const std::vector<double> &primals) {
        serial_ = global_iSQOIterate_serial++;
        primal_values_.assign(&primals[0], &primals[0] + num_primal_);
    }
    void assign_dual_eq(const std::vector<double> &dual_eq) {
        serial_ = global_iSQOIterate_serial++;
        dual_eq_values_.assign(&dual_eq[0], &dual_eq[0] + num_dual_eq_);
    }
    void assign_dual_ieq(const std::vector<double> &dual_ieq) {
        serial_ = global_iSQOIterate_serial++;
        dual_ieq_values_.assign(&dual_ieq[0], &dual_ieq[0] + num_dual_ieq_);
    }
    
    size_t get_serial() const { return serial_; }
// private:
protected:
	size_t num_primal_;
	size_t num_dual_eq_;
	size_t num_dual_ieq_;
	
	std::vector<double> primal_values_;
	std::vector<double> dual_eq_values_;
	std::vector<double> dual_ieq_values_;
	double penalty_parameter_;
    
    size_t serial_;
};

//! \brief print an iterate
//!
//! here, we print out the underlying vectors associated with the iterate iter, to the output stream 'os'.
inline std::ostream& operator<< (std::ostream& os, const iSQOIterate& iter) {
    return iter.print(os);
}

#endif
