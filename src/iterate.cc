
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "step.hh"
#include "iterate.hh"


iSQOIterate::iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq) : 
		num_primal_(number_primal), 
		num_dual_eq_(number_dual_eq),
		num_dual_ieq_(number_dual_ieq),
		primal_values_(number_primal),
		dual_eq_values_(num_dual_eq_),
		dual_ieq_values_(num_dual_ieq_),
		penalty_parameter_(1e-1)
{
	// penalty_parameter_ = 1e-1;
}
iSQOIterate::iSQOIterate(const iSQOIterate& s):
	num_primal_(s.num_primal_), 
	num_dual_eq_(s.num_dual_eq_),
	num_dual_ieq_(s.num_dual_ieq_),
	primal_values_(s.num_primal_),
	dual_eq_values_(s.num_dual_eq_),
	dual_ieq_values_(s.num_dual_ieq_),
	penalty_parameter_(s.penalty_parameter_)
		// bug found here: forgot to initialize the penalty parameter, which fucked up the linesearch.
		{
	// cerr << "copy being made!" << endl;
	
}
double iSQOIterate::x_norm() const {
	double norm=0.0;
	double norm_sqrt=0.0;
	for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		norm += primal_values_[primal_index]*primal_values_[primal_index];
	}
	norm_sqrt = sqrt(norm);
	return norm_sqrt;
}
std::ostream &iSQOIterate::print(std::ostream &os) {
	os << "prim: [";
	for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		if (primal_index != 0) os << ", ";
		os << primal_values_[primal_index];
	}
	os << "];";
	os << "dual eq: [";
	for (size_t dual_eq_index=0; dual_eq_index<num_dual_eq_; ++dual_eq_index) {
		if (dual_eq_index != 0) os << ", ";
		os << dual_eq_values_[dual_eq_index];
	}
	os << "];";
	os << "dual ieq: [";
	for (size_t dual_ieq_index=0; dual_ieq_index<num_dual_ieq_; ++dual_ieq_index) {
		if (dual_ieq_index != 0) os << ", ";
		os << dual_ieq_values_[dual_ieq_index];
		
	}
	os << "];" << std::endl;
}
void iSQOIterate::update(const iSQOIterate &iterate, double alpha, const iSQOStep& step) {
	assert(step.num_primal_ == num_primal_);
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index) {
		primal_values_[primal_index] = iterate.primal_values_[primal_index] + alpha*step.primal_values_[primal_index];
	}
	// for (size_t dual_eq_index=0; dual_eq_index < iterate.num_dual_eq_; ++dual_eq_index) {
	// 	dual_eq_values_[dual_eq_index] = step.dual_eq_values_[dual_eq_index];
	// }
	// for (size_t dual_ieq_index=0; dual_ieq_index < iterate.num_dual_ieq_; ++dual_ieq_index) {
	// 	dual_ieq_values_[dual_ieq_index] = step.dual_ieq_values_[dual_ieq_index];
	// }
}
void iSQOIterate::update_dual(const iSQOStep& step) {
	assert(step.num_dual_eq_ == num_dual_eq_);
	assert(step.num_dual_ieq_ == num_dual_ieq_);
	for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq_; ++dual_eq_index) {
		dual_eq_values_[dual_eq_index] = step.dual_eq_values_[dual_eq_index];
	}
	for (size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq_; ++dual_ieq_index) {
		dual_ieq_values_[dual_ieq_index] = step.dual_ieq_values_[dual_ieq_index];
	}
}
