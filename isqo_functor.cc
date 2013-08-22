// overly-ambitious-iSQO
//  - Travis C. Johnson (traviscj@traviscj.com)
//  - August 2013

// basic philosophy is:
// - keep main loop as clean and as close to algo in the paper as humanly possible, at the expense of anything else.
// - don't try to fit functions (e.g. linear reduction) into classes for iterates
//   - INSTEAD: write functions that stand on their own (but functors, so they can have some persistent state)
// - don't try to fit multiple types of iterates into one class
//   - INSTEAD: separate the iterates from steps from subproblems from hessian shifting from...
// - don't try to fit more functionality than absolutely necessary into classes.
//   - INSTEAD: have easy-to-debug classes with a single simple point: handling hessian shifts, for instance.
// - don't try to pass in the kitchen sink simply to be able to evaluate some basic things.
//   - INSTEAD: Pass the NLP class, which returns f(x), c(x), \bar{c}(x), take it from there.
// - DO NOT use short (1-2 character) or greek symbols as variable names.
//   - INSTEAD: spell out exactly what you are iterating over. Primal indices like "i"? Try primal_index. Linesearch step size called "alpha"? Try linesearch_step_size.
//   - Exception: External class interfaces are a bit messier.
//   - Gotcha: A bit more consistent naming scheme would improve things even more.
// - whenever possible, return ONE object by value. 
//   - INSTEAD: let C++11x worry about RVO. 
//   - Exception: iSQOStep class... but maybe it doesn't need to be?
// - whenever possible, pass required objects by const reference.
// - whenever possible, use vector<double> instead of any other vectors.
// - for unfinished features, use assert to ensure they are not used.

// future feature list:
// - DONE !!! variable bounds !!!
// - DONE !!! more robust shifting strategy and/or Quasi-Newton. (Quasi-Newton might still be an interesting experiment)
// - DONE... rewrite some of the matrix computations into the matrix class. (But still need to check for stragglers...)
// - constraint/objective scaling (should fix HS99... )
// - other iSQO features (algo II-features.)
// - iQP interface! -- but how to do fallback/etc?
// - integrate DJ/CK/AW/et al's work on sparse qpOASES.
// - logbook-style reporting/logging, with sql.
// - use qpoases matrices instead of custom matrix code
// - sparse matrices from ampl/to qpOASES.
// - BQPD interface?
// - read parameters from files, store in some structure.
// - second order correction (might fix HS067...)
// - use multiple qp objects for getting subproblem solutions.
// - const all applicable functions.
// - 

// justification for weirder choices:
// - I use "functors" for some of the paper-defined functions, to capture program state. They work by overloading operator(), which is a big weird, but they are essentially functions with state. One key thing: State is still minimized in them--most are just wrappers around basic NLP functions. With this I was trying to accomplish clear input-output relations, instead of inventing functions that do incremental changes to a huge block of state. One downside to this, as written, is that many things get re-evaluated(This could be minimized by some caching of )

// future ideas:
// - would like to extend the FunctionWithNLPState class to cache the contents of the arguments, and return the last value if they haven't changed. I think this could be somewhat straightforward, but maybe I'm underestimating. This is another major "point" to doing 
// 

// helpful hints:
// - OSX: for alloc bugs, run with: export DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib. On other systems, valgrind.
// - currently, here, the reported KKT is the two-norm... in the matlab implementation, it is the inf norm, scaled by a bunch of things.

// remaining potential bugs:
// hs067.out:Failure - Did not converge - tiny step length... maybe bad step from QP solver? -- might be best fixed by a second-order correction.
// hs099.out:Failure - Did not converge - this *should* be fixed with problem rescaling - should add problem rescaling.
// other HS failures:
// hs087.out:Failure - Did not converge - was also iter limit in MATLAB implmentation.
// hs098.out:Failure - Did not converge - subproblem failure in MATLAB implementation.
// hs106.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs109.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs114.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs117.out:Failure - Did not converge - iteration limit in MATLAB implementation.
// hs99exp.out:Failure - Did not converge - iteration limit in MATLAB implementation.

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>

#include <qpOASES.hpp>

#include "asl.h"

using namespace std;

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector< T >& vec) {
    os << "[";
	for (size_t vector_index = 0; vector_index < vec.size(); ++vector_index) {
		if (vector_index != 0) os << ", ";
		os << vec[vector_index];
	}
    os << " ]";
    return os;
}

double bracket_plus(double val) {
	return max(val, 0.0);
}
double bracket_minus(double val) {
	return max(-val, 0.0);
}
string ordinal(int n) {
	if (n%10==1)
		return string("st");
	if (n%10==2)
		return string("nd");
	if (n%10==3)
		return string("rd");
	else
		return string("th");
}

int isqostep_serial = 0;
class iSQOStep {
public:
	// #define X(s) (std::cout << s << ": " << isqostep_serial << "\n")
	iSQOStep(int number_primal, int number_dual_eq, int number_dual_ieq) : 
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
	iSQOStep(const iSQOStep& other) : serial(isqostep_serial++),
		num_primal_(other.num_primal_),
		num_dual_eq_(other.num_dual_eq_),
		num_dual_ieq_(other.num_dual_ieq_),
		primal_values_(num_primal_),
		dual_eq_values_(num_dual_eq_),
		dual_ieq_values_(num_dual_ieq_)
		{
	    // X("iSQOStep::iSQOStep(const iSQOStep&)");
		// cerr << "copy being made!" << endl;
		num_primal_ = other.num_primal_;
		
		
		for (size_t primal_index=0; primal_index < num_primal_; ++primal_index)
			primal_values_[primal_index] = other.primal_values_[primal_index];
		for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq_; ++dual_eq_index)
			dual_eq_values_[dual_eq_index] = other.dual_eq_values_[dual_eq_index];
		for (size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq_; ++dual_ieq_index)
			dual_ieq_values_[dual_ieq_index] = other.dual_ieq_values_[dual_ieq_index];
		
		status_ = other.status_;
	}
	// ~iSQOStep() {
	//     X("iSQOStep::~iSQOStep()");
	// 	cout << "destructing an isqo step..." << endl;
	// }
	// 
	const iSQOStep *operator=(const iSQOStep& other) {
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
	
	double x_dot_product(vector<double> other_vector) {
		assert(num_primal_ == other_vector.size());
		double dot_product = 0.0;
		for (size_t primal_index=0; primal_index < num_primal_; ++primal_index)
			dot_product += primal_values_[primal_index] * other_vector[primal_index];
		return dot_product;
	}
	double x_norm() const {
		double norm=0.0;
		double norm_sqrt=0.0;
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			norm += primal_values_[primal_index]*primal_values_[primal_index];
		}
		norm_sqrt = sqrt(norm);
		return norm_sqrt;
	}
	
	void print() const {
		cout << "prim: [";
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			if (primal_index != 0) cout << ", ";
			cout << primal_values_[primal_index];
		}
		cout << "];" << endl;
		cout << "dual eq: [";
		for (size_t dual_eq_index=0; dual_eq_index<num_dual_eq_; ++dual_eq_index) {
			if (dual_eq_index != 0) cout << ", ";
			cout << dual_eq_values_[dual_eq_index];
		}
		cout << "];" << endl;
		cout << "dual ieq: [";
		for (size_t dual_ieq_index=0; dual_ieq_index<num_dual_ieq_; ++dual_ieq_index) {
			if (dual_ieq_index != 0) cout << ", ";
			cout << dual_ieq_values_[dual_ieq_index];
			
		}
		cout << "];" << endl;
	}
	
	void convex_combination(const iSQOStep &penalty_step, const iSQOStep &feasibility_step, double combination_step_contribution_from_penalty_step) {
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			primal_values_[primal_index] = combination_step_contribution_from_penalty_step*penalty_step.primal_values_[primal_index] + 
													  (1.0-combination_step_contribution_from_penalty_step)*feasibility_step.primal_values_[primal_index];
		}
	}
	
	void set_primal(const iSQOStep &step) {
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			primal_values_[primal_index] = step.primal_values_[primal_index];
		}
	}
// private:
// protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;

	vector<double> primal_values_;
	vector<double> dual_eq_values_;
	vector<double> dual_ieq_values_;
	
	int status_;
	
	const int serial;
private:
	// const iSQOStep *operator=(const iSQOStep& );
	// iSQOStep(const iSQOStep& s);
};
class iSQOIterate {
public:
	iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq) : 
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
	iSQOIterate(const iSQOIterate& s):
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
	double x_norm() const {
		double norm=0.0;
		double norm_sqrt=0.0;
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			norm += primal_values_[primal_index]*primal_values_[primal_index];
		}
		norm_sqrt = sqrt(norm);
		return norm_sqrt;
	}
	void print() {
		cout << "prim: [";
		for (size_t primal_index=0; primal_index<num_primal_; ++primal_index) {
			if (primal_index != 0) cout << ", ";
			cout << primal_values_[primal_index];
		}
		cout << "];";
		cout << "dual eq: [";
		for (size_t dual_eq_index=0; dual_eq_index<num_dual_eq_; ++dual_eq_index) {
			if (dual_eq_index != 0) cout << ", ";
			cout << dual_eq_values_[dual_eq_index];
		}
		cout << "];";
		cout << "dual ieq: [";
		for (size_t dual_ieq_index=0; dual_ieq_index<num_dual_ieq_; ++dual_ieq_index) {
			if (dual_ieq_index != 0) cout << ", ";
			cout << dual_ieq_values_[dual_ieq_index];
			
		}
		cout << "];" << endl;
	}
	void update(const iSQOIterate &iterate, double alpha, const iSQOStep& step) {
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
	void update_dual(const iSQOStep& step) {
		assert(step.num_dual_eq_ == num_dual_eq_);
		assert(step.num_dual_ieq_ == num_dual_ieq_);
		for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq_; ++dual_eq_index) {
			dual_eq_values_[dual_eq_index] = step.dual_eq_values_[dual_eq_index];
		}
		for (size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq_; ++dual_ieq_index) {
			dual_ieq_values_[dual_ieq_index] = step.dual_ieq_values_[dual_ieq_index];
		}
	}
// private:
// protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;
	
	vector<double> primal_values_;
	vector<double> dual_eq_values_;
	vector<double> dual_ieq_values_;
	double penalty_parameter_;
};
class matrix {
public:
	matrix(size_t rows, size_t columns) : rows_(rows), columns_(columns), data_(rows*columns) {
		// cout << "initializing a matrix..." << endl;
	}
	void set(size_t r, size_t c, double val) {
		data_[columns_*r + c] = val;
	}
	double get(size_t r, size_t c) const {
		return data_[columns_*r + c];
	}
	// double
	void print() const {
		for (size_t r=0; r<rows_; r++) {
			for (size_t c=0; c<columns_; c++) {
				if (c!=0) cout << ", ";
				cout << data_[columns_*r + c];
			}
			if (r+1 != rows_) cout << "; ";
		}
	}
	vector<double> multiply(const vector<double> &x) const {
		assert(x.size() == columns_);
		vector<double> retval(rows_);
		for (size_t current_row=0; current_row<rows_; ++current_row) {
			for (size_t current_col=0; current_col<columns_; ++current_col) {
				retval[current_row] += data_[columns_*current_row + current_col]*x[current_col];
			}
		}
		return retval;
	}
	vector<double> multiply_transpose(const vector<double> &y) const {
		assert(y.size() == rows_);
		// if (y.size() != rows_) throw 20;
		vector<double> retval(columns_);
		for (size_t current_row=0; current_row<rows_; ++current_row) {
			for (size_t current_col=0; current_col<columns_; ++current_col) {
				retval[current_col] += data_[columns_*current_row + current_col]*y[current_row];
			}
		}
		return retval;
	}
	size_t rows_, columns_;
	vector<double> data_;
private:
protected:
	
};
inline std::ostream& operator << (std::ostream& os, const matrix& m) {
    os << "[";
    // for (std::vector<double>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    // {
        // os << " " << *ii;
    // }
	m.print();
    os << " ]";
    return os;
}

class Nlp {
	// this class implements 
public:
	Nlp(int num_primal, int num_dual_eq, int num_dual_ieq) : num_primal_(num_primal), num_dual_eq_(num_dual_eq), num_dual_ieq_(num_dual_ieq) {
		// cout << "-- initializing Nlp" << endl; 
	}
	
	// starting point:
	virtual iSQOIterate initial() = 0;
	
	// zeroth order NLP quantities:
	virtual double objective(const iSQOIterate &iterate) = 0;
	virtual vector<double> constraints_equality(const iSQOIterate &iterate) = 0;
	virtual vector<double> constraints_inequality(const iSQOIterate &iterate) = 0;
	
	// first order NLP quantities:
	virtual vector<double> objective_gradient(const iSQOIterate &iterate) = 0;
	virtual matrix constraints_equality_jacobian(const iSQOIterate &iterate) = 0;
	virtual matrix constraints_inequality_jacobian(const iSQOIterate &iterate) = 0;
	
	// second order NLP quantities:
	virtual matrix lagrangian_hessian(const iSQOIterate &iterate) = 0;
	
	int num_primal() { return num_primal_; }
	int num_dual() {return num_dual_eq_+num_dual_ieq_; }
	int num_dual_eq() { return num_dual_eq_; }
	int num_dual_ieq() { return num_dual_ieq_; }
protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;
};

// TODO: Then, implement it using sparse AMPL instead...
class Hs014 : public Nlp {
public:
	Hs014() : Nlp(2,1,1) {
	}
	
	iSQOIterate initial() {
		iSQOIterate tmp(2,1,1);
		tmp.primal_values_[0] = 2.0;
		tmp.primal_values_[1] = 1.0;
		tmp.dual_eq_values_[0] = -1.0;
		tmp.dual_ieq_values_[0] = 1.0;
		return tmp;
	}
	double objective(const iSQOIterate &iterate) {
		return (iterate.primal_values_[0] - 2)*(iterate.primal_values_[0] - 2) + (iterate.primal_values_[1]-1)*(iterate.primal_values_[1]-1);
	}
	vector<double> constraints_equality(const iSQOIterate &iterate) {
		vector<double> retval(iterate.num_dual_eq_);
		retval[0] = iterate.primal_values_[0] - 2*iterate.primal_values_[1] + 1;
		return retval;
	}
	vector<double> constraints_inequality(const iSQOIterate &iterate) {
		vector<double> retval(iterate.num_dual_ieq_);
		retval[0] = .25*iterate.primal_values_[0]*iterate.primal_values_[0] + iterate.primal_values_[1]*iterate.primal_values_[1] - 1;
		return retval;
	}
	
	vector<double> objective_gradient(const iSQOIterate &iterate) {
		vector<double> ret_gradient(iterate.num_primal_);
		ret_gradient[0] = 2*(iterate.primal_values_[0] - 2);
		ret_gradient[1] = 2*(iterate.primal_values_[1] - 1);
		return ret_gradient;
	}
	matrix constraints_equality_jacobian(const iSQOIterate &iterate){
		matrix ret_jacobian(iterate.num_dual_ieq_, iterate.num_primal_);
		ret_jacobian.set(0,0, 1.0);
		ret_jacobian.set(0,1,-2.0);
		return ret_jacobian;
	}
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate){
		matrix ret_jacobian(iterate.num_dual_ieq_, iterate.num_primal_);
		ret_jacobian.set(0,0,0.5*iterate.primal_values_[0]);
		ret_jacobian.set(0,1,2.0*iterate.primal_values_[1]);
		return ret_jacobian;
	}
	
	matrix lagrangian_hessian(const iSQOIterate &iterate) {
		matrix hess(iterate.num_primal_, iterate.num_primal_);
		if (iterate.penalty_parameter_ > 0) {
			hess.set(0,0,iterate.penalty_parameter_*2.0 + 0.5*iterate.dual_ieq_values_[0]); //   /iterate.penalty_parameter_
			hess.set(1,1,iterate.penalty_parameter_*2.0 + 2.0*iterate.dual_ieq_values_[0]);	//	/iterate.penalty_parameter_
		} else {
			hess.set(0,0,0.5*iterate.dual_ieq_values_[0]);
			hess.set(1,1,2.0*iterate.dual_ieq_values_[1]);
		}
		return hess;
	}
protected:
};

class AmplNlp : public Nlp {
public:
	// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
	AmplNlp(string stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
		if (PRINT_) cout << "Constructing an AmplNlp" << endl;
		ASL *asl;
		asl = ASL_alloc(ASL_read_pfgh);
		nerror_ = new int;
		*nerror_ = 0;
		// char * punt = stub_str.data();
		// if (PRINT_) 
			cout << "Filename: " << stub_str << endl;
		nl_ = jac0dim(&stub_str[0], (stub_str).length());
		if (PRINT_) 
			cout << "ran jac0dim: " << nl_ << endl;
			
		num_primal_ = n_var;
		X0_.resize(num_primal());
		X0 = &X0_[0];
	
		int status = pfgh_read(nl_, ASL_return_read_err | ASL_findgroups);
	
		if (PRINT_) cout << "X0: " << X0[0] << ", " << X0[1] << endl;
		
		for (size_t ampl_constraint_index=0; ampl_constraint_index<n_con; ++ampl_constraint_index) {
			if (PRINT_) cout << "On constraint ampl_constraint_index=" << ampl_constraint_index << ": ";
			bool equality_con = LUrhs[2*ampl_constraint_index]==LUrhs[2*ampl_constraint_index+1] && LUrhs[2*ampl_constraint_index]!=INFINITY;
		
			if (equality_con) {
				if (PRINT_) cout << "equality";
				if (PRINT_) cout << "(push eq)";
				equality_constraints_.push_back(ampl_constraint_index);
			} else if (LUrhs[2*ampl_constraint_index] == -INFINITY && LUrhs[2*ampl_constraint_index+1] > -INFINITY) {
				if (PRINT_) cout << "upper bound";
				if (PRINT_) cout << "(push ineq)";
				inequality_constraints_upper_.push_back(ampl_constraint_index);
			} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY) {
				if (PRINT_) cout << "double-sided";
				if (PRINT_) cout << "(push ineq)";
				if (PRINT_) cout << "(push -ineq)";
				inequality_constraints_upper_.push_back(ampl_constraint_index);
				inequality_constraints_lower_.push_back(ampl_constraint_index);
			} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY) {
				if (PRINT_) cout << "lower bound only";
				if (PRINT_) cout << "(push -ineq)";
				inequality_constraints_lower_.push_back(ampl_constraint_index);
			} else {
				cout << "unsupported";
			}
			if (PRINT_) cout << endl;
		}
		
		
		
		vector<double> x(num_primal());
		// double *x = new double[num_primal()];
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
			x[primal_index] = X0_[primal_index];
		
		bool should_abort = false;
		for (size_t ampl_variable_index =0 ; ampl_variable_index < n_var; ++ampl_variable_index) {
			if (PRINT_) 
				cout << "ampl_variable_index=" << ampl_variable_index << ": l=" << LUv[2*ampl_variable_index] << " <= x_ampl_variable_index= " << x[ampl_variable_index] << " <= u=" << LUv[2*ampl_variable_index+1] << endl;
			bool lower_finite = (LUv[2*ampl_variable_index] != -INFINITY);
			bool upper_finite = (LUv[2*ampl_variable_index+1] != INFINITY);
			bool both_finite = lower_finite && upper_finite;
			bool equal = (LUv[2*ampl_variable_index] == LUv[2*ampl_variable_index+1]);
			if (both_finite && equal) {
				variable_equality_.push_back(ampl_variable_index);
			} else if (both_finite) {
				variable_bound_lower_.push_back(ampl_variable_index);
				variable_bound_upper_.push_back(ampl_variable_index);
			} else if (lower_finite) {
				variable_bound_lower_.push_back(ampl_variable_index);
			} else if (upper_finite) {
				variable_bound_upper_.push_back(ampl_variable_index);
			}			
		}
		
		vector<double> con(n_con);
		conval(&x[0], &con[0], nerror_);
		for (size_t ampl_constraint_index =0 ; ampl_constraint_index < n_con; ++ampl_constraint_index) {
			if (PRINT_) cout << "ampl_variable_index=" << ampl_constraint_index << ": l=" << LUrhs[2*ampl_constraint_index] << " <= c(x_k)= " << con[ampl_constraint_index] << " <= u=" << LUrhs[2*ampl_constraint_index+1] << endl;
		}
		
		num_dual_eq_ = equality_constraints_.size();
		num_dual_ieq_ = inequality_constraints_lower_.size() + inequality_constraints_upper_.size() + variable_bound_lower_.size() + variable_bound_upper_.size();
		asl_ = asl;
	}
	
	~AmplNlp() {
		ASL *asl = asl_;
		
		delete nerror_;
		ASL_free(&asl);
		fclose(nl_);
	}
	
	iSQOIterate initial() {
		iSQOIterate initial(num_primal(), num_dual_eq(), num_dual_ieq());
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			initial.primal_values_[primal_index] = X0_[primal_index];
			if (PRINT_) cout << "initial_" << primal_index << ": " << X0_[primal_index] <<endl;
		}
		vector<double> eq_values = constraints_equality(initial);
		for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq(); ++dual_eq_index) {
			if (eq_values[dual_eq_index] < 0) {
				initial.dual_eq_values_[dual_eq_index] = -1.0;
			} else if (eq_values[dual_eq_index] > 0) {
				initial.dual_eq_values_[dual_eq_index] = +1.0;
			} else {
				initial.dual_eq_values_[dual_eq_index] = 0.0;
			}
		}
	
		vector<double> ieq_values = constraints_inequality(initial);
		for (size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq(); ++dual_ieq_index) {
			if (ieq_values[dual_ieq_index] < 0) {
				initial.dual_ieq_values_[dual_ieq_index] = 0.0;
			} else if (ieq_values[dual_ieq_index] > 0) {
				initial.dual_ieq_values_[dual_ieq_index] =+1.0;
			} else {
				initial.dual_ieq_values_[dual_ieq_index] = 0.0;
			}
		}
		return initial;
	}
	// zeroth order NLP quantities:
	double objective(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> x(num_primal());
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		real obj = objval(0, &x[0], nerror_);
		if (PRINT_) cout << "objective value(nerror = " << *nerror_ << "): " << obj << endl;
		return obj;
	}
	
	vector<double> constraints_equality(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> x(num_primal());
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		// vector<double> con(num_dual());
		vector<double> con(n_con);
		conval(&x[0], &con[0], nerror_);
		if (PRINT_) cout << "equality_constraints:" << endl;
		vector<double> equality_constraint_evaluation(equality_constraints_.size());
		for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
			equality_constraint_evaluation[isqo_eq_constraint_index] = con[equality_constraints_[isqo_eq_constraint_index]] - LUrhs[2*equality_constraints_[isqo_eq_constraint_index]];
		}
		if (PRINT_) cout << "eq eval: " << equality_constraint_evaluation;
		return equality_constraint_evaluation;
	}
	vector<double> constraints_inequality(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> x(num_primal());
		
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		vector<double> con(n_con);
		
		conval(&x[0], &con[0], nerror_);
		if (PRINT_) cout << "inequality_constraints:" << endl;
		
		vector<double> inequality_constraint_evaluation(num_dual_ieq());
		
		int isqo_ineq_constraint_index = 0;
		for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
			if (PRINT_) cout << "isqo ineq " << isqo_ineq_constraint_index << "; " << inequality_constraints_lower_[isqo_ineq_lower_constraint_index] << endl;
			size_t ampl_constraint_index = inequality_constraints_lower_[isqo_ineq_lower_constraint_index];
			size_t ampl_constraint_value_lower = 2*ampl_constraint_index;
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = -con[ampl_constraint_index] + LUrhs[ampl_constraint_value_lower];
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
			
			size_t ampl_constraint_index = inequality_constraints_upper_[isqo_ineq_upper_constraint_index];
			size_t	ampl_constraint_value_upper = 2*ampl_constraint_index+1;
			if (PRINT_) cout << "isqo ineq " << isqo_ineq_upper_constraint_index << " is AMPL constraint " << ampl_constraint_index << endl;
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = con[ampl_constraint_index] - LUrhs[ampl_constraint_value_upper];
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]];
			// cout 
			// 	<< "variable lower bound " << isqo_ineq_lower_variable_index 
			// 	<< "; ampl bound: " <<  variable_bound_lower_[isqo_ineq_lower_variable_index]
			// 	<< "; with value: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]]
			// 	<< "; @x_k: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]]
			// 	<< endl;
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
			// cout << "variable upper bound " << isqo_ineq_upper_variable_index << "" << endl;
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = -LUv[2*variable_bound_lower_[isqo_ineq_upper_variable_index]+1] + x[variable_bound_lower_[isqo_ineq_upper_variable_index]];
			++isqo_ineq_constraint_index;
		}
		assert(isqo_ineq_constraint_index == num_dual_ieq());
		
		
		if (PRINT_) cout <<"ineq eval" << inequality_constraint_evaluation;
		return inequality_constraint_evaluation;
	}
	
	// first order NLP quantities:
	vector<double> objective_gradient(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> return_gradient(num_primal());
		
		vector<double> x(num_primal());
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		objgrd(0, &x[0], &return_gradient[0], nerror_);
		if (PRINT_) cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[0] << endl;
		if (PRINT_) cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[1] << endl;
		
		return return_gradient;
	}
	matrix constraints_equality_jacobian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> x(num_primal());
		
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		matrix equality_constraint_jacobian(equality_constraints_.size(), n_var);
		if (PRINT_) cout << "equal: " << equality_constraints_ << endl;
		vector<double> G(num_primal());
		for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
			congrd(equality_constraints_[isqo_eq_constraint_index], &x[0], &G[0], nerror_);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				equality_constraint_jacobian.set(isqo_eq_constraint_index, primal_index, G[primal_index]);
			}
			
		}
		if (PRINT_) cout << "equality jacobian: [" << endl;
		if (PRINT_) equality_constraint_jacobian.print();
		if (PRINT_) cout << "]" << endl;

		return equality_constraint_jacobian;
	}
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> x(num_primal());
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index){
			x[primal_index] = iterate.primal_values_[primal_index];
		}
		
		matrix inequality_constraint_jacobian(num_dual_ieq(), n_var);
		size_t isqo_ineq_constraint_index=0;
		if (PRINT_) cout << "lower: " << inequality_constraints_lower_ << endl;
		if (PRINT_) cout << "upper: " << inequality_constraints_upper_ << endl;
		vector<double> G(num_primal());
		for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
			congrd(inequality_constraints_lower_[isqo_ineq_lower_constraint_index], &x[0], &G[0], nerror_);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, -G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
			congrd(inequality_constraints_upper_[isqo_ineq_upper_constraint_index], &x[0], &G[0], nerror_);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
			inequality_constraint_jacobian.set(isqo_ineq_constraint_index, variable_bound_lower_[isqo_ineq_lower_variable_index], -1.0);
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
			inequality_constraint_jacobian.set(isqo_ineq_constraint_index, variable_bound_upper_[isqo_ineq_upper_variable_index], +1.0);
			++isqo_ineq_constraint_index;
		}
		if (PRINT_) 
		{
			cout << "======" << endl;
			cout << inequality_constraint_jacobian << endl;
			cout << "======" << endl;
		}
				
		return inequality_constraint_jacobian;
	}
	
	// second order NLP quantities:
	matrix lagrangian_hessian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		vector<double> H(num_primal() * num_primal());
		vector<double> OW(1);
		OW[0] = iterate.penalty_parameter_;
		vector<double> Y(n_con);
		
		// funnel dual_{,i}eq_values_ into Y to pass to AMPL.
		for (size_t dual_eq_index=0; dual_eq_index < iterate.dual_eq_values_.size(); ++dual_eq_index) {
			Y[equality_constraints_[dual_eq_index]] = iterate.dual_eq_values_[dual_eq_index];// / iterate.penalty_parameter_;
		}
		
		size_t dual_ineq_value_index = 0;
		for (size_t dual_ineq_lower_index=0; dual_ineq_lower_index < inequality_constraints_lower_.size(); ++dual_ineq_lower_index) {
			Y[inequality_constraints_lower_[dual_ineq_lower_index]] = -iterate.dual_ieq_values_[dual_ineq_value_index];
			dual_ineq_value_index++;
		}
		for (size_t dual_ineq_upper_index=0; dual_ineq_upper_index < inequality_constraints_upper_.size(); ++dual_ineq_upper_index) {
			Y[inequality_constraints_upper_[dual_ineq_upper_index]] = iterate.dual_ieq_values_[dual_ineq_value_index];
			dual_ineq_value_index++;
		}
		if (PRINT_) cout << "Y[0] = " << Y[0] << endl;
		if (PRINT_) cout << "Y[1] = " << Y[1] << endl;
		if (PRINT_) cout << "Y[2] = " << Y[2] << endl;
		if (PRINT_) cout << "Y[3] = " << Y[3] << endl;		
		
		// In this call:
		//	 - H : OUTPUT : vector holding entries of full hessian.
		//	 - n_var : INPUT : # of variables or something about the stride, maybe.
		//	 - 0 : INPUT : index of the desired objective.
		//	 - OW : INPUT : multipliers for objective function ("objective weights")
		//	 - Y : INPUT : lagrange multipliers for constraints.
		fullhes(&H[0], n_var, 0, &OW[0], &Y[0]);
		if (iterate.penalty_parameter_ == 0) {
			vector<double> H_f_only(num_primal()*num_primal());
			vector<double> OW_f_only(1);
			OW_f_only[0] = 1.0;
			vector<double> Y_f_only(n_con);
			for (size_t dual_index=0; dual_index < n_con; ++dual_index) Y[dual_index] = 0.0;
			fullhes(&H_f_only[0], n_var, 0, &OW_f_only[0], &Y_f_only[0]);
			
			for (size_t row_index=0; row_index < n_var; ++row_index)
				for (size_t col_index=0; col_index < n_var; ++col_index) 
					H[row_index*n_var + col_index] -= H_f_only[row_index*n_var + col_index];
		}
		
		matrix return_hessian(n_var, n_var);
		for (size_t row_index=0; row_index < n_var; ++row_index) {
			for (size_t col_index=0; col_index < n_var; ++col_index) {
				return_hessian.set(row_index, col_index, H[row_index*n_var + col_index]);
			}
		}
		return return_hessian;
	}
	
	
	size_t num_lower_ieq() {
		return inequality_constraints_lower_.size();
	}
	size_t num_upper_ieq() {
		return inequality_constraints_upper_.size();
	}
private:
protected:
	bool PRINT_;
	ASL *asl_;
	FILE *nl_;
	
	vector<double> X0_;
	fint *nerror_;
	vector<size_t> equality_constraints_;
	vector<size_t> inequality_constraints_lower_;
	vector<size_t> inequality_constraints_upper_;
	
	vector<size_t> variable_equality_;
	vector<size_t> variable_bound_lower_;
	vector<size_t> variable_bound_upper_;
};

class FunctionWithNLPState {
public:
	FunctionWithNLPState(Nlp &nlp) {
		nlp_ = &nlp;
	}
protected:
	Nlp *nlp_;
};

class ConstraintViolationFunction : public FunctionWithNLPState{
public:
	ConstraintViolationFunction(Nlp &nlp) : FunctionWithNLPState(nlp){
		// cout << "-- Initializing l1 violation function." << endl;
	}
	double operator()(const iSQOIterate &iterate) const {
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		
		return two_vectors(con_values_eq, con_values_ieq);
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
		bool PRINT = false;

		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);

		if (PRINT) cout << "constraintviolationfunction e 0: " << con_values_eq << endl;
		// J(x)*d
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		if (PRINT) jac_ce.print();
		for (size_t row=0; row<iterate.num_dual_eq_; ++row) {
			for (size_t col=0; col<iterate.num_primal_; ++col) {
				con_values_eq[row] += jac_ce.get(row,col)*step.primal_values_[col];
			}
		}
		if (PRINT) cout << "constraintviolationfunction e 1: " << con_values_eq << endl;
		
		// \bar{J}(x)*d
		if (PRINT) cout << "constraintviolationfunction i 0: " << con_values_ieq << endl;
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		// if (PRINT) jac_ci.print();
		for (size_t row=0; row<iterate.num_dual_ieq_; ++row) {
			for (size_t col=0; col<iterate.num_primal_; ++col) {
				con_values_ieq[row] += jac_ci.get(row,col)*step.primal_values_[col];
			}
		}
		if (PRINT) cout << "constraintviolationfunction i 1: " << con_values_ieq << endl;
		
		return two_vectors(con_values_eq, con_values_ieq);
	}
	
	double two_vectors(const vector<double> &eq_items, const vector<double> &ieq_items) const {
		double eq_violation = 0.0;
		for (size_t dual_eq_index=0; dual_eq_index<eq_items.size(); ++dual_eq_index) {
			eq_violation += abs(eq_items[dual_eq_index]);
		}
		double ieq_violation = 0.0;
		for (size_t dual_ieq_index=0; dual_ieq_index<ieq_items.size(); ++dual_ieq_index) {
			// if (ieq_items[i] > 0){
			ieq_violation += bracket_plus(ieq_items[dual_ieq_index]);
			// }
		}
		return eq_violation + ieq_violation;
	}
protected:
};

class PenaltyFunction : public FunctionWithNLPState{
public:
	PenaltyFunction(Nlp &nlp) : 
			FunctionWithNLPState(nlp), 
			constraint_violation_func_(nlp)
			{
		// cout << "-- Initializing penalty function." << endl;
	}
	double operator()(const iSQOIterate &iterate) const {
		// cout << "--- Calling the PenaltyFunction Functor..." << endl;
		double f = nlp_->objective(iterate);
		// cout << "--- objective: " << f << endl;
		return iterate.penalty_parameter_*f + constraint_violation_func_(iterate);
	}
protected:
	ConstraintViolationFunction constraint_violation_func_;
	// double penalty_parameter_;
};

class iSQOQuadraticSubproblem : public FunctionWithNLPState{
public:
	iSQOQuadraticSubproblem(Nlp &nlp, const iSQOIterate &iterate) :
				FunctionWithNLPState(nlp),
				num_qp_variables_(nlp.num_primal() + 2*nlp.num_dual()), 
				num_qp_constraints_(nlp.num_dual()),
				num_nlp_variables_(nlp.num_primal()), 
				num_nlp_constraints_eq_(nlp.num_dual_eq()), 
				num_nlp_constraints_ieq_(nlp.num_dual_ieq()),
				hessian_(num_qp_variables_, num_qp_variables_), 
				unshifted_hessian_diagonal_(num_qp_variables_), 
				first_shift_(false),
				jacobian_(num_qp_constraints_, num_qp_variables_), 
				gradient_(num_qp_variables_),
				lower_bound_(num_qp_variables_), 
				upper_bound_(num_qp_variables_), 
				jacobian_lower_bound_(num_qp_constraints_), 
				jacobian_upper_bound_(num_qp_constraints_),
				nlp_hessian_(num_nlp_variables_, num_nlp_variables_),
				nlp_eq_jacobian_(num_nlp_constraints_eq_, num_nlp_variables_),
				nlp_ieq_jacobian_(num_nlp_constraints_ieq_, num_nlp_variables_),
				nlp_objective_gradient_(num_nlp_variables_)
	 {
		nlp_objective_gradient_ = nlp_->objective_gradient(iterate);
		nlp_hessian_ = nlp_->lagrangian_hessian(iterate);
		nlp_eq_jacobian_ = nlp_->constraints_equality_jacobian(iterate);
		nlp_ieq_jacobian_ = nlp_->constraints_inequality_jacobian(iterate);
	
		for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(eq_constraint_index,variables, nlp_eq_jacobian_.get(eq_constraint_index,variables));
			}
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+eq_constraint_index, -1.0);
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+iterate.num_dual_eq_+eq_constraint_index, +1.0);
		}
				
		for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,variables, nlp_ieq_jacobian_.get(ieq_constraint_index,variables));
			}
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+ieq_constraint_index, -1.0);
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+ieq_constraint_index, +1.0);
		}		
		for (size_t r=0; r<iterate.num_primal_; ++r) {
			for (size_t c=0; c<iterate.num_primal_; ++c) {
				hessian_.set(r,c, nlp_hessian_.get(r,c));
			}
		}
				
		// NLP gradient copied over
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index)
			gradient_[primal_index] = iterate.penalty_parameter_*nlp_objective_gradient_[primal_index];
		// equality slacks both have positive:
		for (size_t dual_eq_index=0; dual_eq_index<iterate.num_dual_eq_; ++dual_eq_index) {
			gradient_[iterate.num_primal_+dual_eq_index] = 1.0;
			gradient_[iterate.num_primal_+iterate.num_dual_eq_+dual_eq_index] = 1.0;
		}
		// inequality slacks only penalize the positive parts:
		for (size_t dual_ieq_index=0; dual_ieq_index<iterate.num_dual_ieq_; ++dual_ieq_index) {
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+dual_ieq_index] = 1.0;
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+dual_ieq_index] = 0.0;
		}
		
		vector<double> con_values_eq=nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);
		
		for (size_t eq_index=0; eq_index<iterate.num_dual_eq_; ++eq_index) {
			jacobian_lower_bound_[eq_index] = -con_values_eq[eq_index];
			jacobian_upper_bound_[eq_index] = -con_values_eq[eq_index];
		}
		for (size_t ieq_index=0; ieq_index<iterate.num_dual_ieq_; ++ieq_index) {
			jacobian_lower_bound_[iterate.num_dual_eq_+ieq_index] = -INFINITY;
			jacobian_upper_bound_[iterate.num_dual_eq_+ieq_index] = -con_values_ieq[ieq_index];
		}
		for (size_t variable_index=0; variable_index < iterate.num_primal_; ++variable_index) {
			lower_bound_[variable_index] = -1e10;
			upper_bound_[variable_index] = +1e10;
		}
		for (size_t variable_index=0; variable_index < 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_); ++variable_index) {
			lower_bound_[iterate.num_primal_+variable_index] = 0.0;
			upper_bound_[iterate.num_primal_+variable_index] = 1e10;
		}
	}
	
	void inc_regularization(double hessian_shift) {
		bool PRINT = false;
		// make a copy of the diagonal entries, if we haven't already got one...
		if (! first_shift_) {
			if (PRINT) cout << "making the first-time-only copy..." << endl;
			for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
					unshifted_hessian_diagonal_[diagonal_entry] = hessian_.get(diagonal_entry,diagonal_entry);
			}
			first_shift_ = true;
		}
		
		for (size_t diagonal_entry=0; diagonal_entry<num_nlp_variables_; ++diagonal_entry) {
				hessian_.set(diagonal_entry,diagonal_entry, unshifted_hessian_diagonal_[diagonal_entry] + hessian_shift);
		}
		for (size_t diagonal_entry=num_nlp_variables_; diagonal_entry < num_qp_variables_; ++diagonal_entry) {
			hessian_.set(diagonal_entry, diagonal_entry, hessian_shift);
		}
	}
	
	void print() const {
		cout << endl;
		cout << "H = [";
		hessian_.print();
		cout << "];" << endl;
		cout << "A = [";
		jacobian_.print();
		cout << "];" << endl;
		
		cout << "g = " << gradient_ << "';" << endl;
		cout << "lb = " << lower_bound_ << "';" << endl;
		cout << "ub = " << upper_bound_ << "';" << endl;
		cout << "lbA = " << jacobian_lower_bound_ << "';" << endl;
		cout << "ubA = " << jacobian_upper_bound_ << "';" << endl;
	}
	int num_qp_variables_, num_qp_constraints_;
	int num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_;
	matrix hessian_;
	bool first_shift_;
	vector<double> unshifted_hessian_diagonal_;
	matrix jacobian_;
	vector<double> gradient_;
	vector<double> lower_bound_;
	vector<double> upper_bound_;
	vector<double> jacobian_lower_bound_;
	vector<double> jacobian_upper_bound_;
	matrix nlp_hessian_;
	matrix nlp_eq_jacobian_;
	matrix nlp_ieq_jacobian_;
	vector<double> nlp_objective_gradient_;
private:
protected:
};

class ResidualFunction : public FunctionWithNLPState {
public:
	ResidualFunction(Nlp &nlp) : FunctionWithNLPState(nlp) {
		// cout << "-- Initializing nlp residual function." << endl;
	}
	double resid_helper(const iSQOIterate &iterate, vector<double> stationarity, vector<double> constraint_eq_values, vector<double> constraint_ieq_values, vector<double> constraint_eq_dual_values, vector<double> constraint_ieq_dual_values) const {
		vector<double> rho(stationarity.size() + 2*constraint_eq_values.size() + 2*constraint_ieq_values.size());
		
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
				
		size_t rho_index = 0;

		double eq_signflip = 1;
		double ieq_signflip = 1;
		// cout << "sta pre: " << stationarity << endl << endl << endl;
		vector<double> eq_jacobian_trans_times_mults = jac_ce.multiply_transpose(constraint_eq_dual_values);
		vector<double> ieq_jacobian_trans_times_mults = jac_ci.multiply_transpose(constraint_ieq_dual_values);
		
		for (size_t stationarity_index=0; stationarity_index < stationarity.size(); ++stationarity_index) {
			rho[rho_index] = stationarity[stationarity_index] + eq_jacobian_trans_times_mults[stationarity_index] + ieq_jacobian_trans_times_mults[stationarity_index];
			++rho_index;
		}
		for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {
			rho[rho_index] = min(bracket_plus(constraint_eq_values[constraint_eq_index]), 1.0 - eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
			
			rho[rho_index] = min(bracket_minus(constraint_eq_values[constraint_eq_index]), 1.0 + eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
		}
		
		for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
			rho[rho_index] = min(bracket_plus(constraint_ieq_values[constraint_ieq_index]), 1.0 - ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
			++rho_index;

			rho[rho_index] = min(bracket_minus(constraint_ieq_values[constraint_ieq_index]), 0.0 + ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
			++rho_index;
			
		}
		// cout << "rho(" << iterate.penalty_parameter_ << "): " << rho << endl;
		double total = 0.0;
		for (size_t rho_index = 0; rho_index < rho.size(); ++rho_index)
			total += rho[rho_index]*rho[rho_index];
		
		return sqrt(total);
	}
	double operator()(const iSQOIterate &iterate) const {
		bool PRINT=false;
		// cout << "OPERATOR FOR RESID FUNC: " << iterate.penalty_parameter_ << endl;
		
		vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f:
			stationarity[stationarity_index] *= iterate.penalty_parameter_;
		}
		
		vector<double> eq_con_values = nlp_->constraints_equality(iterate);
		vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
		
		if (PRINT) {
			cout << endl << "mu: " << iterate.penalty_parameter_ << endl;
			cout << "sta (pre): " << nlp_->objective_gradient(iterate) << endl;
			cout << "sta: " << stationarity << endl;
			cout << "constraints  eq: " << eq_con_values << endl;
			cout << "constraints ieq: " << ieq_con_values << endl;
		
			cout << "jac eq: " << nlp_->constraints_equality_jacobian(iterate) << endl;
		
			cout << "       dual  eq: " << iterate.dual_eq_values_ << endl;
			cout << "       dual ieq: " << iterate.dual_ieq_values_ << endl;
		}
		return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, iterate.dual_eq_values_, iterate.dual_ieq_values_);
	}
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// first entry of rho(x,y,bary,mu):
		bool PRINT=true;
		vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		matrix hessian = nlp_->lagrangian_hessian(iterate);
		
		vector<double> hessian_times_step = hessian.multiply(step.primal_values_);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f + H*d:
			stationarity[stationarity_index] = iterate.penalty_parameter_*stationarity[stationarity_index] + hessian_times_step[stationarity_index];
			
		}
		
		matrix jacobian_eq = nlp_->constraints_equality_jacobian(iterate);
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		// cout << "about to do j eq" << endl;
		vector<double> linear_step_eq = jacobian_eq.multiply(step.primal_values_);
		
		// cout << "about to do j ieq" << endl;
		matrix jacobian_ieq = nlp_->constraints_inequality_jacobian(iterate);
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		vector<double> linear_step_ieq = jacobian_ieq.multiply(step.primal_values_);

		for (size_t eq_con_index=0; eq_con_index< nlp_->num_dual_eq(); ++eq_con_index) {
			con_values_eq[eq_con_index] += linear_step_eq[eq_con_index];
		}
		for (size_t ieq_con_index=0; ieq_con_index< nlp_->num_dual_ieq(); ++ieq_con_index) {
			con_values_ieq[ieq_con_index] += linear_step_ieq[ieq_con_index];
		}		
		return resid_helper(iterate, stationarity, con_values_eq, con_values_ieq, step.dual_eq_values_, step.dual_ieq_values_);
	}
protected:
};

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual()), first_(true)  {}
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem) {
		// qpOASES::SQProblem example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual());
		/* Setting up QProblem object. */
		int nWSR = 2000;
		qpOASES::returnValue ret;
		// cout << endl;
		// cout << "hessian: " << subproblem.hessian_ <<endl;
		// cout << "jacobian: " << subproblem.jacobian_ <<endl;
		// cout << "gradient: " << subproblem.gradient_ << endl;
		// cout << "lb: : " << subproblem.lower_bound_ << endl;
		// cout << "ub: : " << subproblem.upper_bound_ << endl;
		// cout << "lbA: " << subproblem.jacobian_lower_bound_ << endl;
		// cout << "ubA: " << subproblem.jacobian_upper_bound_ << endl;
		
		qpOASES::Options opt;
		// opt.setToReliable();
		// opt.enableRegularisation = qpOASES::BooleanType(true);

		opt.terminationTolerance = 1e-6;
		example_.setOptions(opt);
		example_.setPrintLevel(qpOASES::PL_NONE);
		// cout.flush();
		// cerr.flush();
		if (first_) {
			// cout << "subproblem.hessian_.data_[0]: " << subproblem.hessian_ << endl;
			// cout << "subproblem: " << subproblem.gradient_ << endl;
			// cout << "subproblem: " << subproblem.jacobian_ << endl;
			// cout << "subproblem: " << subproblem.lower_bound_ << endl;
			// cout << "subproblem: " << subproblem.upper_bound_ << endl;
			// cout << "subproblem: " << subproblem.jacobian_lower_bound_ << endl;
			// cout << "subproblem: " << subproblem.jacobian_upper_bound_ << endl;
			// subproblem.print();
			ret = example_.init( &subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR,
								 0);
		} else{
			
			// subproblem.print();
			ret = example_.hotstart(&subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR,
								 0);
		}
		cerr.flush();
		cout.flush();

		size_t return_status = example_.getStatus();
		// getGlobalMessageHandler()->listAllMessages();
		if( ret != qpOASES::SUCCESSFUL_RETURN ){
			// subproblem.print();
	        // printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
			first_ = true;
			example_.reset();
		} else {
			if (first_) first_ = false;
		}
		iSQOStep step(subproblem.num_nlp_variables_,subproblem.num_nlp_constraints_eq_,subproblem.num_nlp_constraints_ieq_);
		// cout << "address of step: " << &step << endl;
		
		// full_primal needs to hold primal variables and all slack variables
		vector<double> primal_and_slack_values(nlp_->num_primal() + 2*nlp_->num_dual());
		example_.getPrimalSolution( &primal_and_slack_values[0] );
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
			step.primal_values_[primal_index] = primal_and_slack_values[primal_index];
		
		// cout << "primal: " << step.primal_values_ << endl;
		// 	- every QP variable has a multiplier:
		//		-- iterate.num_primal_: primal variables in NLP
		// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
		//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
		//	- every constraint has a multiplier:
		//		-- iterate.num_dual_eq_
		//		-- iterate.num_dual_ieq_
		// so in total: 
		vector<double> yOpt(subproblem.num_qp_variables_ + subproblem.num_qp_constraints_);
		example_.getDualSolution( &yOpt[0] );				
		
		step.status_ = ret;
		// to get eq con multipliers, read past NLP & slack variables:
		for (size_t dual_eq_index=0; dual_eq_index<subproblem.num_nlp_constraints_eq_; ++dual_eq_index) {
			step.dual_eq_values_[dual_eq_index] = -yOpt[subproblem.num_qp_variables_+dual_eq_index];
		}
		// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
		for (size_t dual_ieq_index=0; dual_ieq_index<subproblem.num_nlp_constraints_ieq_; ++dual_ieq_index) {
			step.dual_ieq_values_[dual_ieq_index] = -yOpt[subproblem.num_qp_variables_+subproblem.num_nlp_constraints_eq_+dual_ieq_index];
		}
		// cout << "dual eq: " << step.dual_eq_values_ << endl;
		// cout << "dual ieq: " << step.dual_ieq_values_ << endl;
		
		bool PRINT=false;
		if (PRINT) {
			step.print();
		}
		return step;
	}
private:
protected:
	qpOASES::SQProblem example_;
	bool first_;
};

class LinearModelFunction : public FunctionWithNLPState {
public:
	LinearModelFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
		
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
		double f = nlp_->objective(iterate);
		vector<double> gradient = nlp_->objective_gradient(iterate);
		// dotproduct(gradient, )
		double dot_product=0.0;
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index) {
			dot_product += gradient[primal_index]*step.primal_values_[primal_index];
		}
		return iterate.penalty_parameter_*(f + dot_product) + constraint_violation_func_(iterate, step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

class LinearReductionFunction : public FunctionWithNLPState {
public:
	LinearReductionFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
		
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
		bool PRINT=false;
		// double f = nlp_->objective(iterate);
		vector<double> gradient = nlp_->objective_gradient(iterate);
		// dotproduct(gradient, )
		double dot_product=0.0;
		for (size_t primal_index=0; primal_index<iterate.num_primal_; ++primal_index) {
			dot_product += gradient[primal_index]*step.primal_values_[primal_index];
		}
		
		if (PRINT) cout << "linear decrease: " << endl;
		if (PRINT) cout << " - " << -iterate.penalty_parameter_*dot_product << endl;
		if (PRINT) cout << " - " << constraint_violation_func_(iterate) << endl;
		if (PRINT) cout << " - " << constraint_violation_func_(iterate,step) << endl;
		if (PRINT) cout << " - " << -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step) << endl;
		return -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

class LineSearchFunction : public FunctionWithNLPState {
public:
	LineSearchFunction(Nlp &nlp) : FunctionWithNLPState(nlp), penalty_func_(nlp), linear_decrease_func_(nlp) {
		// cout << "initializing a linesearcher!" << endl;
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
		// bugs found here: it_and_step pulled from nlp->initial, which meant it didn't have the right penalty parameter, which caused the linesearch to fail.
		iSQOIterate it_and_step(iterate);
		bool PRINT = false;
		double step_size = 1.0;
		double penfunc_start = penalty_func_(iterate);
		double linear_reduction_start = linear_decrease_func_(iterate, step);
		double machine_precision = std::numeric_limits<double>::epsilon();
		
		double eta = 1e-8;
		double linesearch_step_size_reduction = 0.5;
		if (PRINT) cout << endl;
		for (size_t alpha_cuts = 0; alpha_cuts < 20; alpha_cuts++) {
			it_and_step.update(iterate, step_size, step);
			double penfunc_step = penalty_func_(it_and_step);
			
			if (PRINT) {
				cout << "step_size: " << step_size << endl;
				cout << "f: " << nlp_->objective(it_and_step) << endl;
				cout << "c: " << nlp_->constraints_equality(it_and_step) << endl;
				cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << endl;
				cout << "  rhs: " << - eta*step_size*linear_reduction_start  << " + " << 10*machine_precision*max(abs(penfunc_step), 1.0) << " = " << (-eta*step_size*linear_reduction_start+10*machine_precision*max(abs(penfunc_step), 1.0)) << endl;
				cout << "mu = " << it_and_step.penalty_parameter_ << endl;
			}
			if (penfunc_step - penfunc_start <= - eta*step_size*linear_reduction_start + 10*machine_precision*max(abs(penfunc_step), 1.0)) {
				break;
			}
			step_size = linesearch_step_size_reduction*step_size;
		}
		return step_size;
	}
private:
protected:
	PenaltyFunction penalty_func_;
	LinearReductionFunction linear_decrease_func_;
};
class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(Nlp &nlp) : FunctionWithNLPState(nlp), solve_qp_(nlp), last_shift_(0.0) {
		
	}
	iSQOStep operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem) {
		bool PRINT=false;
		double shift_w0 = 1e-4;
		double shift_min = 1e-20;
		double shiftkwbarp = 100;
		double shiftkwp = 8;
		double shift_kwm = 1.0/3.0;
		double shiftmax = 1e40;
		
		double current_shift = 0.0;
		
		// try to solve without regularization:
		subproblem.inc_regularization(current_shift);
		iSQOStep return_step = solve_qp_(subproblem);
		
		
		int current_regularization_steps = 0;
		vector<double> current_step_values(nlp_->num_primal() + 2*nlp_->num_dual());
		for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
			current_step_values[primal_index] = return_step.primal_values_[primal_index];
		}
		
		// matrix hessian = nlp_->lagrangian_hessian(iterate);
		vector<double> hessian_step = subproblem.hessian_.multiply(current_step_values);
		if (PRINT) cout << "hessian step: " << subproblem.hessian_ << endl;
		double total = 0.0;
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
			total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
		}
		if (PRINT) cout << endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << endl;
		if (PRINT) cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << endl;
		// vector<double> hessian_step = subproblem.hessian_.multiply()
		while ((return_step.status_ != 0) || (return_step.x_norm() > 1e9) || (total < (1e-8*return_step.x_norm()*return_step.x_norm()))) {
			if (current_regularization_steps == 0) {
				if (last_shift_ == 0.0) {
					current_shift = shift_w0;
				} else {
					current_shift = max(shift_min, shift_kwm * last_shift_);
				}
			} else {
				if (last_shift_ == 0.0) {
					current_shift = shiftkwbarp * current_shift;
				} else {
					current_shift = shiftkwp*current_shift;
				}
			}
			subproblem.inc_regularization(current_shift);
			return_step = solve_qp_(subproblem);
			
			for (size_t primal_index = 0; primal_index < nlp_->num_primal(); ++primal_index) {
				current_step_values[primal_index] = return_step.primal_values_[primal_index];
			}
			if (PRINT) cout << "hessian: " << subproblem.hessian_ << endl;
			if (PRINT) cout << "return step: " << return_step.primal_values_ << endl;
			hessian_step = subproblem.hessian_.multiply(current_step_values);
			if (PRINT) cout << "hessian step: " << hessian_step << endl;
			total = 0.0;
			for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index) {
				total += hessian_step[primal_index] * return_step.primal_values_[primal_index];
			}
			if (PRINT) cout << "current_shift: " << current_shift << endl;
			if (PRINT) cout << endl << "total: " << total << "; norm: " << 1e-8*return_step.x_norm()*return_step.x_norm() << endl;
			if (PRINT) cout << "required: " << (total < (1e-8*return_step.x_norm()*return_step.x_norm())) << endl;
			
			++current_regularization_steps;
			
			if (current_regularization_steps == 20){
				cout << "FAILURE!" << endl;
				assert(false);
				break;
			}
		}
		last_shift_ = current_shift;
		return return_step;
	}
	
	double get_last_shift() const {
		return last_shift_;
	}
private:
protected:
	SolveQuadraticProgram solve_qp_;
	double last_shift_;
};

class TextOutput : public FunctionWithNLPState {
public:
	TextOutput (Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp),linear_decrease_func_(nlp), pen_func_(nlp), residual_func_(nlp) {
		
	}
	
	void start() {
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-----------------[ FEASIBILITY ]-----------------|");
		printf("-------------------[ PENALTY ]-------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
		printf("%s",output_desc_pre_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_subprob_);
		printf("%s",output_desc_post_);
		printf("-----|");
		printf("----------------------|");
		printf("----------------------|");
		printf("----------------------&");
		printf("-------------------------------------------------|");
		printf("-------------------------------------------------|");
		printf("--------------------------------------------|");
		printf("----------\n");
	}
	void pre(size_t iter, const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate) const {
		// iter, problem, feas iterate, pen iterate, constraintviolation, iterate, penfunc, residual_func
		//  -> params: iter, feas iter, pen iter
		//	-> state: nlp_, constraint_violation_func_, pen_func_, residual_func_
		printf(output_format_pre_, 
				iter, 
				nlp_->objective(penalty_iterate), constraint_violation_func_(penalty_iterate),
				penalty_iterate.penalty_parameter_, pen_func_(penalty_iterate),
				residual_func_(feasibility_iterate), residual_func_(penalty_iterate)
				);
	}
	void subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// shift, step, linear_decrease_func, residual_func
		//  -> params: shift, step
		//	-> state: linear_decrease_func, residual_func_
		printf(output_format_subprob_,
				shift, step.status_, step.x_norm(),
				linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
					);
	}
	void subproblem_skip() {
		printf("         -   -           -          -          - |");
	}
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, string step_type, double step_mix, double alpha) {
		// step, feas_iter, penalty_iter, alpha
		//  -> params: step, feas_iter, penalty_iter, alpha
		printf(output_format_post_,
				step_type.c_str(),
				step_mix, combination_step.x_norm(), linear_decrease_func_(feasibility_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
				alpha
					);
	}
protected:
	// Nlp *nlp_;
	ConstraintViolationFunction constraint_violation_func_;
	LinearReductionFunction linear_decrease_func_;
	PenaltyFunction pen_func_;
	ResidualFunction residual_func_;
	
	static const char output_desc_pre_[];
	static const char output_desc_subprob_[];
	static const char output_desc_post_[];
	static const char output_format_pre_[];
	static const char output_format_subprob_[];
	static const char output_format_post_[];
};
const char TextOutput::output_desc_pre_[] = " it  |       obj     infeas |       pen      merit |   feaskkt     penkkt &";
const char TextOutput::output_desc_subprob_[] = "     shift  msg     ||d||     penred        res  |";
const char TextOutput::output_desc_post_[] = " TT   CvxComb    ||d||   FeasRed     PenRed |    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |%9.2e\n";


bool assert_close(double val1, double val2, double tol) {
	assert(abs(val1) <= abs(val2) + tol*abs(val1));
	assert(abs(val2) <= abs(val1) + tol*abs(val2));
	assert(((val1 > tol) && (val2 > tol)) || ((val1 < -tol) && (val2 < -tol)) || (val1 == 0 && val2 == 0));
	return (abs(val1) <= abs(val2) + tol*abs(val1)) && (abs(val2) <= abs(val1) + tol*abs(val2));
}
int main3() {
	string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs009.nl");
	AmplNlp problem(problem_file);
	iSQOIterate penalty_iterate = problem.initial();
	matrix je = problem.constraints_equality_jacobian(penalty_iterate);
	vector<double> x(penalty_iterate.primal_values_);
	vector<double> y(penalty_iterate.dual_eq_values_);
	cout << "jac: "<< je << endl;
	cout << "x: " << x << endl;
	cout << je.multiply(x) << endl;
	
	cout << "y: " << y << endl;
	cout << je.multiply_transpose(y) << endl;
	
}
int main2() {
	cout << "okay, so..." << endl;
	
	// tests for constraint conversation:
	// equality constrained problems: (to c(x) = 0)
	//	* 0 <= c_1(x) <= 0  ---> A_1d = 0 - c_1(x_k)
	//	* 1 <= c_2(x) <= 1  ---> A_2d = 1 - c_2(x_k)
	// inequality constrained problems: (to barc(x) <= 0)
	//  * double constrained: 	-1 <= c_3(x) <= 5  --> A_3d <= 5 - c_3(x_k) AND -A_3d <= -(-1) - (-c_3(x_k))
	//			c_3(x) <= 5  ==> 
	//			-1 <= c_3(x) ==> -c_3(x) <= -(-1) ==> -A_3d - c_3(x_k) <= -(-1) ==> -A_3d  <= -(-1) - (-c_3(x_k))
	//	* double constrained: 	0 <= c_4(x) <= 5 --> A_4d <= 5 - c_4(x_k) AND -A_4d <= -(0) - (-c_4(x_k))
	//	* double constrained: 	-1 <= c_5(x) <= 0 --> A_5d <= 0 - c_5(x_k) AND -A_5d <= -(-1) - (-c_5(x_k))
	//	* only upper: -INF <= c_6(x) <= 0 --> A_6d <= 0 - c_6(x_k)
	//	* only upper: -INF <= c_7(x) <= 5 --> A_7d <= 5 - c_7(x_k)
	//	* only lower: 0 <= c_8(x) <= INF --> -A_8d <= -(0) - c_8(x_k)
	//	* only lower: 5 <= c_9(x) <= INF --> -A_9d <= -(5) - c_9(x_k)
	
	// 
	// tests for variable conversion:
	// 
	// AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs016.nl");
	AmplNlp ampl_nlp_object("/Users/traviscj/ampl_bound_reformulation.nl");
	iSQOIterate iterate = ampl_nlp_object.initial();
	cout << "initial iterate: " << iterate.primal_values_ << endl;
	// cout << "Zeroth Order:" << endl;
	cout << "===[Objective]==============" << endl;
	double ampl_obj = ampl_nlp_object.objective(iterate);
	assert(ampl_obj == 5.921590000000000e+00);//, "ampl_obj wrong");
	cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << endl;
	
	cout << "===[Objective Gradient]==============" << endl;
	vector<double> ampl_obj_gradient = ampl_nlp_object.objective_gradient(iterate);
	assert(ampl_obj_gradient[0] == 1.0e0);
	assert(ampl_obj_gradient[1] == 1.0e0);
	cout << "  AMPL : " << ampl_obj_gradient << endl;
	
	
	cout << "===[Equality Constraints]==============" << endl;
	vector<double> con_eq_values = ampl_nlp_object.constraints_equality(iterate);
	cout << "  AMPL : " << con_eq_values << endl;
	assert(con_eq_values[0]==-41.924375456199996); // negative because switches var to left-hand-side.
	assert(con_eq_values[1]==141.09951409669998);
	
	cout << "===[Equality Jacobian]==============" << endl;
	matrix eq_jacobian = ampl_nlp_object.constraints_equality_jacobian(iterate);
	cout << "  AMPL : " << eq_jacobian << endl;
	assert(eq_jacobian.get(0,0) == -1.256636000000000e+01);
	assert(eq_jacobian.get(0,1) == -1.668000000000000e+01);
	assert(eq_jacobian.get(1,0) == 4.398226000000000e+01);
	assert(eq_jacobian.get(1,1) == 6.116000000000000e+01);
	
	
	
	
	cout << "===[Inequality Constraints]==============" << endl;
	vector<double> ieq_values = ampl_nlp_object.constraints_inequality(iterate);
	cout << "  AMPL : " << ieq_values << endl;
	double ieq_tol = 1e-12;
	assert(abs(ieq_values[0]) <= (abs(-3.482753668339000e+02) + ieq_tol*abs(ieq_values[0])));
	assert(abs(ieq_values[1]) <= (abs(-6.820391459396999e+02) + ieq_tol*abs(ieq_values[1])));
	assert(abs(ieq_values[2]) <= (abs(3.362753668339000e+02) + ieq_tol*abs(ieq_values[2])));
	assert(abs(ieq_values[3]) <= (abs(7.346270723082998e+02) + ieq_tol*abs(ieq_values[3])));
	
	// -6.820391459396999e+02
	
	cout << "===[Inequality Jacobian]==============" << endl;
	matrix ieq_jacobian = ampl_nlp_object.constraints_inequality_jacobian(iterate);
	cout << "  AMPL : " << ieq_jacobian << endl;
	assert_close(ieq_jacobian.get(0,0), -1.193804200000000e+02, 1e-10); // c2 lower
	assert_close(ieq_jacobian.get(0,1), -1.278800000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(1,0), -2.324776600000000e+02, 1e-10); // c3 lower
	assert_close(ieq_jacobian.get(1,1), -2.279600000000000e+02, 1e-10);
	assert_close(ieq_jacobian.get(2,0), 1.193804200000000e+02, 1e-10);  // c2 upper
	assert_close(ieq_jacobian.get(2,1), 1.278800000000000e+02, 1e-10);  
	assert_close(ieq_jacobian.get(3,0), 2.701767400000000e+02, 1e-10);  // c4 upper
	assert_close(ieq_jacobian.get(3,1), 2.613200000000000e+02, 1e-10);
	
	
	
	
	
	
	// assert(ieq_jacobian.get(3,0) == 1.193804200000000e+02); assert(ieq_jacobian.get(3,1) == 1.278800000000000e+02);  // c2 upper
	
	cout << endl;
	cout << "TESTED FRONTIER!" << endl;
	cout << endl;
	cout << endl;
	
	cout << "===[Lagrangian Hessian]==============" << endl;
	matrix lh_f_only = ampl_nlp_object.lagrangian_hessian(iterate);
	cout << "  AMPL : " << lh_f_only << endl;
	assert_close(lh_f_only.get(0,0), 0.00e0, 1e-10);  // c2 upper
	assert_close(lh_f_only.get(0,1), 0.00e0, 1e-10);  
	assert_close(lh_f_only.get(1,0), 0.00e0, 1e-10);  // c4 upper
	assert_close(lh_f_only.get(1,1), 0.00e0, 1e-10);
	
	iterate.dual_eq_values_[0] = 1.0;
	matrix lh_c0_only = ampl_nlp_object.lagrangian_hessian(iterate);
	cout << "  AMPL : " << lh_c0_only << endl;
	iterate.dual_eq_values_[0] = -1.0;
	matrix lh_c0_only_flip = ampl_nlp_object.lagrangian_hessian(iterate);
	cout << "  AMPL : " << lh_c0_only_flip << endl;
	
	
	return 0;
}
int main1() {
	AmplNlp ampl_nlp_object("/Users/traviscj/optimization/cute_nl_nopresolve/hs014.nl");
	
	
	Hs014 problem;
	iSQOIterate iterate(2,1,1);
	iterate.penalty_parameter_ = 1.0e-1;
	iterate.primal_values_[0] = 2.0;
	iterate.primal_values_[1] = 2.0;
	iterate.dual_eq_values_[0] = -1.0;
	iterate.dual_ieq_values_[0] = 1000;
	// iterate.dual_ieq_values_[0] = ;
	// iterate.dual_ieq_values_[0] = 0;
	// iterate.penalty_parameter_ = 1.0e-1;
	// iterate.primal_values_[0] = .822875656;
	// iterate.primal_values_[1] = .911437828;
	// iterate.dual_eq_values_[0] = 1.846589027861980e+00;
	// iterate.dual_ieq_values_[0] = 1.594493103554523e+00;

	cout << "Zeroth Order:" << endl;
	cout << "===[Objective]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.objective(iterate) << endl;
	cout << "Custom : " << problem.objective(iterate) << endl;
	
	cout << "===[Equality Constraints]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.constraints_equality(iterate) << endl;
	cout << "Custom : " << problem.constraints_equality(iterate) << endl;
	
	cout << "===[Inequality Constraints]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.constraints_inequality(iterate) << endl;
	cout << "Custom : " << problem.constraints_inequality(iterate) << endl;
	
	cout << "First Order:" << endl;
	cout << "===[Objective Gradient]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.objective_gradient(iterate) << endl;
	cout << "Custom : " << problem.objective_gradient(iterate) << endl;
	
	cout << "===[Equality Jacobian]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.constraints_equality_jacobian(iterate) << endl;
	cout << "Custom : " << problem.constraints_equality_jacobian(iterate) << endl;
	
	
	cout << "===[Inequality Jacobian]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.constraints_inequality_jacobian(iterate) << endl;
	cout << "Custom : " << problem.constraints_inequality_jacobian(iterate) << endl;
	
	cout << "Second Order:" << endl;
	cout << "===[Lagrangian Hessian]==============" << endl;
	cout << "  AMPL : " << ampl_nlp_object.lagrangian_hessian(iterate) << endl;
	cout << "Custom : " << problem.lagrangian_hessian(iterate) << endl;
		
	return 0;
}
int main(int argc, char **argv) {
	
	string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/airport.nl");
	if (argc>1) {
		problem_file = string(argv[1]);
	}
	AmplNlp problem(problem_file);
	
	// Utilities for NLP:
	PenaltyFunction penalty_function(problem);
	ResidualFunction kkt_residual(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraint_violation(problem);
	LineSearchFunction line_search(problem);
	LinearReductionFunction linear_reduction(problem);
	HessianShifter hessian_shifting_penalty_qp_solve(problem);
	HessianShifter hessian_shifting_feasibility_qp_solve(problem);
	// Other
	TextOutput text_output(problem);
	
	// Paper-defined Parameter Values
	double linear_decrease_threshold = 1e-1; // epsilon in the paper.
	double linear_reduction_threshold_for_penalty_reduction = 1e-2; // beta in the paper.
	double mostly_feasibility_step_threshold = 1e-3; // tau in the paper.
	double penalty_parameter_reduction_factor = 0.2; // delta in paper.
	double termination_threshold = 1e-6;
	double termination_penalty_threshold = 1e-8;
	double hessian_min_convexity = 1e-8;	// theta in the paper.
	// Code-defined parameter values:
	double convex_combination_search_decrease = .9;
	bool PRINT=false;
	size_t maximum_iterations = 1000;
	size_t max_num_comb_reductions = 100;
	double machine_precision = std::numeric_limits<double>::epsilon();
	
	text_output.start();
	
	//////////////////////////
	// ALGORITHM A // STEP 1
	//////////////////////////
	iSQOIterate penalty_iterate = problem.initial();
	penalty_iterate.penalty_parameter_ = 1.0e-1;

	iSQOIterate feasibility_iterate = problem.initial();
	feasibility_iterate.penalty_parameter_ = 0.0;
	size_t iter=-1;
	for (iter = 0; iter < maximum_iterations; iter ++ ) {
		iSQOStep combination_step(problem.num_primal(), problem.num_dual_eq(), problem.num_dual_ieq());
		
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		
		//////////////////////////
		// ALGORITHM A // STEP 2
		//////////////////////////
		if ((kkt_residual(penalty_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) < termination_threshold)) {
			cout << endl << "Termination 2a - optimality" << endl;
			break;
		}
		if ((kkt_residual(feasibility_iterate) < termination_threshold) && (constraint_violation(penalty_iterate) > termination_threshold) && (penalty_iterate.penalty_parameter_ < termination_penalty_threshold)) {
			cout << endl << "Termination 2b - infeasible problem!" << endl;
			break;
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 3
		//////////////////////////
		// Penalty problem is set up AND SOLVED, 
		iSQOQuadraticSubproblem penalty_subproblem(problem, penalty_iterate);
		iSQOStep penalty_step = hessian_shifting_penalty_qp_solve(penalty_iterate, penalty_subproblem);
		iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
		
		// per-iteration variable setup.
		string steptype = "4a";
		double combination_step_contribution_from_penalty_step = -1.0;
		
		if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*constraint_violation(penalty_iterate) + 10*machine_precision*linear_reduction(penalty_iterate,penalty_step)) {
			//////////////////////////
			// ALGORITHM A // STEP 3a
			combination_step_contribution_from_penalty_step = 1.0;
			text_output.subproblem_skip();

			combination_step.set_primal(penalty_step);
		} else {
			
			//////////////////////////
			// ALGORITHM A // STEP 4
			//////////////////////////
			// Feasibility problem is set up and solved.
			iSQOQuadraticSubproblem feasibility_subproblem(problem, feasibility_iterate);
			feasibility_step = hessian_shifting_feasibility_qp_solve(feasibility_iterate, feasibility_subproblem);
			
			if (linear_reduction(penalty_iterate,penalty_step) >= linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step)) {
				//////////////////////////
				// ALGORITHM A // STEP 4a
				combination_step_contribution_from_penalty_step = 1.0;
				steptype = "5a";
			} else {
				//////////////////////////
				// ALGORITHM A // STEP 4b
				combination_step_contribution_from_penalty_step = 1.0;
				
				int num_comb_reductions = 0;
				combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
				while (! (linear_reduction(feasibility_iterate,combination_step)  >=linear_decrease_threshold*linear_reduction(feasibility_iterate,feasibility_step) + 10*machine_precision*linear_reduction(feasibility_iterate,combination_step)) ) {
					combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
					
					combination_step_contribution_from_penalty_step = convex_combination_search_decrease*combination_step_contribution_from_penalty_step;
					++num_comb_reductions;
					assert(! (num_comb_reductions >= max_num_comb_reductions));
				}
				steptype = "5b";
				
				// beta: 
				double linear_decrease_in_penalty_combination = linear_reduction(penalty_iterate,combination_step);
				double linear_decrease_in_feasibility_combination = linear_reduction(feasibility_iterate,combination_step);
				if (linear_decrease_in_penalty_combination >= linear_reduction_threshold_for_penalty_reduction*linear_decrease_in_feasibility_combination + 10*machine_precision*linear_decrease_in_penalty_combination) {
					if (combination_step_contribution_from_penalty_step >= mostly_feasibility_step_threshold) {
						// no-op: keep the current penalty parameter.
					} else {
						penalty_iterate.penalty_parameter_ = penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_;
					}
				} else {
					vector<double> gradient = problem.objective_gradient(penalty_iterate);
					double potential_new_penalty_parameter_numerator = (1-linear_reduction_threshold_for_penalty_reduction)*linear_decrease_in_feasibility_combination;
					double potential_new_penalty_parameter_denominator = combination_step.x_dot_product(gradient) + hessian_min_convexity*pow(combination_step.x_norm(),2);
					penalty_iterate.penalty_parameter_ = min(penalty_parameter_reduction_factor*penalty_iterate.penalty_parameter_, potential_new_penalty_parameter_numerator / potential_new_penalty_parameter_denominator);
				}
			}

			// combination gets the convex combination of penalty and feasibility:
			combination_step.convex_combination(penalty_step, feasibility_step, combination_step_contribution_from_penalty_step);
			
			text_output.subproblem(hessian_shifting_feasibility_qp_solve.get_last_shift(), feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		double step_size = line_search(penalty_iterate, combination_step);
		
		text_output.subproblem(hessian_shifting_penalty_qp_solve.get_last_shift(), penalty_iterate, penalty_subproblem, penalty_step);
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, combination_step_contribution_from_penalty_step, step_size);
		
		//////////////////////////
		// ALGORITHM A // STEP 6
		//////////////////////////
		// update penalty iterate & dual values:
		penalty_iterate.update(penalty_iterate, step_size, combination_step);
		penalty_iterate.update_dual(penalty_step);
		// update feasibility iterate & dual values:
		feasibility_iterate.update(penalty_iterate, step_size, combination_step);
		feasibility_iterate.update_dual(feasibility_step);
	}
	
	if (iter == maximum_iterations) {
		cout << "Failure - Did not converge" << endl;
	}
	cout << penalty_iterate.primal_values_ << endl;
	return 0;
}


