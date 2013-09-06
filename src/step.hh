
#ifndef __GUARD_STEP_HH
#define __GUARD_STEP_HH

#include <iostream>
#include <vector>

class iSQOStep {
public:
	iSQOStep(int number_primal, int number_dual_eq, int number_dual_ieq);
	iSQOStep(const iSQOStep& other);
	const iSQOStep *operator=(const iSQOStep& other);
	
	double x_dot_product(std::vector<double> other_vector);
	double x_norm() const;
	std::ostream &print(std::ostream&) const;
	void convex_combination(const iSQOStep &penalty_step, const iSQOStep &feasibility_step, double combination_step_contribution_from_penalty_step);
	void set_primal(const iSQOStep &step);
    
// private:
// protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;

	std::vector<double> primal_values_;
	std::vector<double> dual_eq_values_;
	std::vector<double> dual_ieq_values_;
	
	int status_;
	
	const int serial;
private:
};

std::ostream& operator<<(std::ostream& os, const iSQOStep &step);

#endif
