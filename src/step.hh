// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#ifndef __GUARD_STEP_HH
#define __GUARD_STEP_HH

#include <iostream>
#include <vector>

//! \brief class for storing iSQO steps (primal values and one set of dual values)
class iSQOStep {
public:
	iSQOStep(int number_primal, int number_dual_eq, int number_dual_ieq, int status, int pivots);
	iSQOStep(const iSQOStep& other);
	const iSQOStep *operator=(const iSQOStep& other);
	
	double x_dot_product(std::vector<double> other_vector);
	double x_norm() const;
	std::ostream &print(std::ostream&) const;
	void convex_combination(const iSQOStep &penalty_step, const iSQOStep &feasibility_step, double combination_step_contribution_from_penalty_step);
	void set_primal(const iSQOStep &step);
    
    void set_primal_value(int index, double value)  { primal_values_[index] = value; }
    void set_dual_eq_value(int index, double value)  { dual_eq_values_[index] = value; }
    void set_dual_ieq_value(int index, double value)  { dual_ieq_values_[index] = value; }
    double get_primal_value(int index) const { return primal_values_[index]; }
    double get_dual_eq_value(int index) const { return dual_eq_values_[index]; }
    double get_dual_ieq_value(int index) const { return dual_ieq_values_[index]; }
    const std::vector<double> &get_primal_values() const { return primal_values_; }
    const std::vector<double> &get_dual_eq_values() const { return dual_eq_values_; }
    const std::vector<double> &get_dual_ieq_values() const { return dual_ieq_values_; }
    

    int get_pivots() const { return pivots_; }    
    void set_pivots(int num_pivots) { pivots_ = num_pivots; }
    int get_status() const { return status_; }
    void set_status(int status) { status_ = status; }
    
    int num_primal() const { return num_primal_; }
    int num_dual_eq() const { return num_dual_eq_; }
    int num_dual_ieq() const { return num_dual_ieq_; }
    
// private:
protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;

	std::vector<double> primal_values_;
	std::vector<double> dual_eq_values_;
	std::vector<double> dual_ieq_values_;

	int pivots_;	
	int status_;
	
	const int serial_;
private:
};

std::ostream& operator<<(std::ostream& os, const iSQOStep &step);

#endif
