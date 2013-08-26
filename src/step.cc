
#include <iostream>

#include <cassert>
#include <cmath>

#include "step.hh"

int isqostep_serial = 0;

// #define X(s) (std::cout << s << ": " << isqostep_serial << "\n")
iSQOStep::iSQOStep(int number_primal, int number_dual_eq, int number_dual_ieq) : 
			num_primal_(number_primal), 
			num_dual_eq_(number_dual_eq),
			num_dual_ieq_(number_dual_ieq),
			primal_values_(number_primal),
			dual_eq_values_(num_dual_eq_),
			dual_ieq_values_(num_dual_ieq_),
			status_(-42),
			serial(isqostep_serial++)
	{
        // X("iSQOStep::iSQOStep()");
	}
iSQOStep::iSQOStep(const iSQOStep& other) : serial(isqostep_serial++),
	num_primal_(other.num_primal_),
	num_dual_eq_(other.num_dual_eq_),
	num_dual_ieq_(other.num_dual_ieq_),
	primal_values_(num_primal_),
	dual_eq_values_(num_dual_eq_),
	dual_ieq_values_(num_dual_ieq_)
	{
    // X("iSQOStep::iSQOStep(const iSQOStep&)");
	// cerr << "copy being made!" << std::endl;
	num_primal_ = other.num_primal_;
	
	
	for (std::size_t primal_index=0; primal_index < num_primal_; ++primal_index)
		primal_values_[primal_index] = other.primal_values_[primal_index];
	for (std::size_t dual_eq_index=0; dual_eq_index < num_dual_eq_; ++dual_eq_index)
		dual_eq_values_[dual_eq_index] = other.dual_eq_values_[dual_eq_index];
	for (std::size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq_; ++dual_ieq_index)
		dual_ieq_values_[dual_ieq_index] = other.dual_ieq_values_[dual_ieq_index];
	
	status_ = other.status_;
}
// ~iSQOStep() {
//     X("iSQOStep::~iSQOStep()");
// 	cout << "destructing an isqo step..." << std::endl;
// }
// 
const iSQOStep *iSQOStep::operator=(const iSQOStep& other) {
	iSQOStep tmp(other);
	
	std::swap(num_primal_, tmp.num_primal_);
	std::swap(num_dual_eq_, tmp.num_dual_eq_);
	std::swap(num_dual_ieq_, tmp.num_dual_ieq_);
	
	std::swap(primal_values_, tmp.primal_values_);
	std::swap(dual_eq_values_, tmp.dual_eq_values_);
	std::swap(dual_ieq_values_, tmp.dual_ieq_values_);
	
	std::swap(status_, tmp.status_);
	return this;
}

double iSQOStep::x_dot_product(std::vector<double> other_vector) {
	assert(num_primal_ == other_vector.size());
	double dot_product = 0.0;
	for (std::size_t primal_index=0; primal_index < num_primal_; ++primal_index)
		dot_product += primal_values_[primal_index] * other_vector[primal_index];
	return dot_product;
}
double iSQOStep::x_norm() const {
	double norm=0.0;
	double norm_sqrt=0.0;
	for (std::size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		norm += primal_values_[primal_index]*primal_values_[primal_index];
	}
	norm_sqrt = sqrt(norm);
	return norm_sqrt;
}

std::ostream& iSQOStep::print(std::ostream& os) const {
	os << "prim: [";
	for (std::size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		if (primal_index != 0) os << ", ";
		os << primal_values_[primal_index];
	}
	os << "];" << std::endl;
	os << "dual eq: [";
	for (std::size_t dual_eq_index=0; dual_eq_index<num_dual_eq_; ++dual_eq_index) {
		if (dual_eq_index != 0) os << ", ";
		os << dual_eq_values_[dual_eq_index];
	}
	os << "];" << std::endl;
	os << "dual ieq: [";
	for (std::size_t dual_ieq_index=0; dual_ieq_index<num_dual_ieq_; ++dual_ieq_index) {
		if (dual_ieq_index != 0) os << ", ";
		os << dual_ieq_values_[dual_ieq_index];
		
	}
	os << "];" << std::endl;
}

void iSQOStep::convex_combination(const iSQOStep &penalty_step, const iSQOStep &feasibility_step, double combination_step_contribution_from_penalty_step) {
	for (std::size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		primal_values_[primal_index] = combination_step_contribution_from_penalty_step*penalty_step.primal_values_[primal_index] + 
												  (1.0-combination_step_contribution_from_penalty_step)*feasibility_step.primal_values_[primal_index];
	}
}

void iSQOStep::set_primal(const iSQOStep &step) {
	for (std::size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
		primal_values_[primal_index] = step.primal_values_[primal_index];
	}
}

std::ostream& operator<<(std::ostream& os, const iSQOStep &step) {
    step.print(os);
}
