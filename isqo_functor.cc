// overly-ambitious-iSQO
//  - Travis C. Johnson (traviscj@traviscj.com)
//  - August 2013

#include <iostream>
#include <vector>
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
			// penalty_parameter_ = 1e-1;
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
		
	double x_norm() const {
		double norm=0.0;
		double norm_sqrt=0.0;
		for (size_t i=0; i<num_primal_; ++i) {
			// cout << "i=" << i << ": " << primal_values_[i] << endl;
			norm += primal_values_[i]*primal_values_[i];
		}
		norm_sqrt = sqrt(norm);
		return norm_sqrt;
	}
	
	void print() {
		cout << "prim: [";
		for (size_t i=0; i<num_primal_; ++i) {
			// cout << "i=" << i << ": " << primal_values_[i] << endl;
			if (i!=0) cout << ", ";
			cout << primal_values_[i];
		}
		cout << "];" << endl;
		cout << "dual eq: [";
		for (size_t i=0; i<num_dual_eq_; ++i) {
			if (i!=0) cout << ", ";
			cout << dual_eq_values_[i];
		}
		cout << "];" << endl;
		cout << "dual ieq: [";
		for (size_t i=0; i<num_dual_ieq_; ++i) {
			if (i!=0) cout << ", ";
			cout << dual_ieq_values_[i];
			
		}
		cout << "];" << endl;
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
		dual_ieq_values_(s.num_dual_ieq_)
			{
		cerr << "copy being made!" << endl;
		
	}
	double x_norm() const {
		double norm=0.0;
		double norm_sqrt=0.0;
		for (size_t i=0; i<num_primal_; ++i) {
			norm += primal_values_[i]*primal_values_[i];
		}
		norm_sqrt = sqrt(norm);
		return norm_sqrt;
	}
	void print() {
		cout << "prim: [";
		for (size_t i=0; i<num_primal_; ++i) {
			// cout << "i=" << i << ": " << primal_values_[i] << endl;
			if (i!=0) cout << ", ";
			cout << primal_values_[i];
		}
		cout << "];";// << endl;
		cout << "dual eq: [";
		for (size_t i=0; i<num_dual_eq_; ++i) {
			if (i!=0) cout << ", ";
			cout << dual_eq_values_[i];
		}
		cout << "];";// << endl;
		cout << "dual ieq: [";
		for (size_t i=0; i<num_dual_ieq_; ++i) {
			if (i!=0) cout << ", ";
			cout << dual_ieq_values_[i];
			
		}
		cout << "];" << endl;
	}
	void update(const iSQOIterate &iterate, double alpha, const iSQOStep& step) {
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
	matrix(int rows, int columns) : rows_(rows), columns_(columns), data_(rows*columns) {
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
	// vector<double> multiply(const vector<double> x) {
	// 	// for (size_t i=0; i<)
	// 	vector<double> retval(3);
	// 	retval[0] = 1;
	// 	retval[1] = 2;
	// 	retval[2] = 3;
	// 	return retval;
	// }
	int rows_, columns_;
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
		// char * punt = stub_str.data();
		cout << "Filename: " << stub_str << endl;
		nl_ = jac0dim(&stub_str[0], (stub_str).length());
		// if (PRINT_) 
			cout << "ran jac0dim: " << nl_ << endl;
			
		num_primal_ = n_var;
		X0_.resize(num_primal());
		X0 = &X0_[0];
	
		int status = pfgh_read(nl_, ASL_return_read_err | ASL_findgroups);
	
		// for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
	// 		X0_[primal_index] = X0[primal_index];
	// 	}
		
		// if (PRINT_) cout << "ran pfgh_read: " << endl;
		// if (PRINT_) cout << "status: " << status << endl;	
		// 
		// if (PRINT_) cout << "n_var: " << n_var << endl;
		// if (PRINT_) cout << "n_con: " <<n_con << endl;
		// if (PRINT_) if (A_vals != NULL) {
		// 	cout << "A_vals " << A_rownos << endl;
		// 	cout << "A_vals " << A_colstarts << endl;
		// 	cout << "A_vals: " << A_vals[0] << endl;
		// 	cout << "A_vals: " << A_vals[1] << endl;
		// 	cout << "A_vals: " << A_vals[2] << endl;
		// 	cout << "A_vals: " << A_vals[3] << endl;
		// }
	
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
		
		
		nerror_ = new int;
		vector<double> x(num_primal());
		// double *x = new double[num_primal()];
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
			x[primal_index] = X0_[primal_index];
		for (size_t ampl_variable_index =0 ; ampl_variable_index < n_var; ++ampl_variable_index) {
			if (PRINT_) 
				cout << "i=" << ampl_variable_index << ": l=" << LUv[2*ampl_variable_index] << " <= x_i= " << x[ampl_variable_index] << " <= u=" << LUv[2*ampl_variable_index+1] << endl;
			// TODO add variable bounds here, using similar logic to the above.
			// for now, though.... just punt.
			assert(LUv[2*ampl_variable_index] == -INFINITY);
			assert(LUv[2*ampl_variable_index+1] == INFINITY);
		}
		
		// double *con = new double[n_con];
		vector<double> con(n_con);
		conval(&x[0], &con[0], nerror_);
		for (size_t ampl_constraint_index =0 ; ampl_constraint_index < n_con; ++ampl_constraint_index) {
			if (PRINT_) cout << "i=" << ampl_constraint_index << ": l=" << LUrhs[2*ampl_constraint_index] << " <= c(x_k)= " << con[ampl_constraint_index] << " <= u=" << LUrhs[2*ampl_constraint_index+1] << endl;
		}
		
		num_dual_eq_ = equality_constraints_.size();
		num_dual_ieq_ = inequality_constraints_lower_.size() + inequality_constraints_upper_.size();
		asl_ = asl;
	}
	
	~AmplNlp() {
		delete nerror_;
		free(asl_); asl_ = NULL;
		// free(nl_); nl_ = NULL;
		fclose(nl_);
	}
	
	iSQOIterate initial() {
		iSQOIterate initial(num_primal(), num_dual_eq(), num_dual_ieq());
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			initial.primal_values_[primal_index] = X0_[primal_index];
			cout << "initial_" << primal_index << ": " << X0_[primal_index] <<endl;
		}
		vector<double> eq_values = constraints_equality(initial);
		for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq(); ++dual_eq_index) {
			if (eq_values[dual_eq_index] < 0) {
				initial.dual_eq_values_[dual_eq_index] = -1.0; // -1 in iSQO-matlab
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
				initial.dual_ieq_values_[dual_ieq_index] = 0.0;	// 0 in isqo-matlab
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
		// delete[] x;
		// delete[] con;
		return equality_constraint_evaluation;
	}
	vector<double> constraints_inequality(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		// double *x = new double[num_primal_];
		vector<double> x(num_primal());
		
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		// double *con = new double[n_con];
		vector<double> con(num_dual());
		
		conval(&x[0], &con[0], nerror_);
		if (PRINT_) cout << "inequality_constraints:" << endl;
		
		vector<double> inequality_constraint_evaluation(num_dual_ieq());
		
		// int ampl_index = 0;
		int isqo_ineq_constraint_index = 0;
		for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
			if (PRINT_) cout << "isqo ineq " << isqo_ineq_constraint_index << "; " << inequality_constraints_lower_[isqo_ineq_lower_constraint_index] << endl;
			int ampl_constraint_index = inequality_constraints_lower_[isqo_ineq_lower_constraint_index];
			int	ampl_constraint_value_lower = 2*ampl_constraint_index;
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = -con[ampl_constraint_index] + LUrhs[ampl_constraint_value_lower];
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
			
			int ampl_constraint_index = inequality_constraints_upper_[isqo_ineq_upper_constraint_index];
			int	ampl_constraint_value_upper = 2*ampl_constraint_index+1;
			if (PRINT_) cout << "isqo ineq " << isqo_ineq_upper_constraint_index << " is AMPL constraint " << ampl_constraint_index << endl;
			inequality_constraint_evaluation[isqo_ineq_constraint_index] = con[ampl_constraint_index] - LUrhs[ampl_constraint_value_upper];
			++isqo_ineq_constraint_index;
		}
		
		if (PRINT_) cout <<"ineq eval" << inequality_constraint_evaluation;
		// delete[] x;
		// delete[] con;
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
		// double *x = new double[num_primal_];
		vector<double> x(num_primal());
		
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		matrix equality_constraint_jacobian(equality_constraints_.size(), n_var);
		if (PRINT_) cout << "equal: " << equality_constraints_ << endl;
		// double *G = new double[num_primal_];
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

		// delete[] G;
		// delete[] x;
		return equality_constraint_jacobian;
	}
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate) {
		// PRINT_=false;
		ASL *asl = asl_;
		// double *x = new double[num_primal_];
		vector<double> x(num_primal());
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index){
			x[primal_index] = iterate.primal_values_[primal_index];
			// cout << "x[" << primal_index << "]: " << x[primal_index] << endl;
		}
		
		matrix inequality_constraint_jacobian(num_dual_ieq(), n_var);
		size_t isqo_ineq_constraint_index=0;
		if (PRINT_) cout << "lower: " << inequality_constraints_lower_ << endl;
		if (PRINT_) cout << "upper: " << inequality_constraints_upper_ << endl;
		// double *G = new double[num_primal_];
		vector<double> G(num_primal());
		for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
			// cout << ""
			congrd(inequality_constraints_lower_[isqo_ineq_lower_constraint_index], &x[0], &G[0], nerror_);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				// cout << "J_{" << isqo_ineq_lower_constraint_index << ", " << primal_index << "} = " << G[primal_index] << endl;
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, -G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
			congrd(inequality_constraints_upper_[isqo_ineq_upper_constraint_index], &x[0], &G[0], nerror_);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				// cout << "J_{" << isqo_ineq_upper_constraint_index << ", " << primal_index << "} = " << G[primal_index] << endl;
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		
		if (PRINT_) 
		{
			cout << "======" << endl;
			cout << inequality_constraint_jacobian << endl;
			cout << "======" << endl;
		}
		
		// delete[] x;
		// delete[] G;
		
		if (PRINT_) cout << "inequality: [" << endl;
		if (PRINT_) inequality_constraint_jacobian.print();
		if (PRINT_) cout << "]" << endl;
		// PRINT_=false;
		return inequality_constraint_jacobian;
	}
	
	// second order NLP quantities:
	matrix lagrangian_hessian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		// double *H = new double[n_var*n_var];
		vector<double> H(num_primal() * num_primal());
		// double *OW = new double[1];
		vector<double> OW(1);
		OW[0] = iterate.penalty_parameter_;
		vector<double> Y(n_con);
		// double *Y = new double[n_con];
		
		// funnel dual_{,i}eq_values_ into Y to pass to AMPL.
		for (size_t dual_eq_index=0; dual_eq_index < iterate.dual_eq_values_.size(); ++dual_eq_index) {
			Y[equality_constraints_[dual_eq_index]] = iterate.dual_eq_values_[dual_eq_index];// / iterate.penalty_parameter_;
		}
		
		size_t dual_ineq_value_index = 0;
		for (size_t dual_ineq_lower_index=0; dual_ineq_lower_index < inequality_constraints_lower_.size(); ++dual_ineq_lower_index) {
			Y[inequality_constraints_lower_[dual_ineq_lower_index]] = -iterate.dual_ieq_values_[dual_ineq_value_index];// / iterate.penalty_parameter_;
			dual_ineq_value_index++;
		}
		for (size_t dual_ineq_upper_index=0; dual_ineq_upper_index < inequality_constraints_upper_.size(); ++dual_ineq_upper_index) {
			Y[inequality_constraints_upper_[dual_ineq_upper_index]] = iterate.dual_ieq_values_[dual_ineq_value_index];// / iterate.penalty_parameter_;
			dual_ineq_value_index++;
		}
		if (PRINT_) cout << "Y[0] = " << Y[0] << endl;
		if (PRINT_) cout << "Y[1] = " << Y[1] << endl;
		if (PRINT_) cout << "Y[2] = " << Y[2] << endl;
		if (PRINT_) cout << "Y[3] = " << Y[3] << endl;
		// cout << "Y[0] = " << Y[0] << endl;
		
		
		// In this call:
		//	 - H : OUTPUT : vector holding entries of full hessian.
		//	 - n_var : INPUT : # of variables or something about the stride, maybe.
		//	 - 0 : INPUT : index of the desired objective.
		//	 - OW : INPUT : multipliers for objective function ("objective weights")
		//	 - Y : INPUT : lagrange multipliers for constraints.
		fullhes(&H[0], n_var, 0, &OW[0], &Y[0]);
		if (iterate.penalty_parameter_ == 0) {
			// cout << "Oops!" << endl;
			// double *H_f_only = new double[num_primal_];
			vector<double> H_f_only(num_primal()*num_primal());
			// double *OW_f_only = new double[num_primal_];
			vector<double> OW_f_only(1);
			OW_f_only[0] = 1.0;
			// double *Y_f_only = new double[num_primal_];
			vector<double> Y_f_only(n_con);
			for (size_t dual_index=0; dual_index < n_con; ++dual_index) Y[dual_index] = 0.0;
			fullhes(&H_f_only[0], n_var, 0, &OW_f_only[0], &Y_f_only[0]);
			
			for (size_t row_index=0; row_index < n_var; ++row_index)
				for (size_t col_index=0; col_index < n_var; ++col_index) 
					H[row_index*n_var + col_index] -= H_f_only[row_index*n_var + col_index];
		}
		
		if (PRINT_) cout << "H(nerror = " << *nerror_ << "): " << H[0] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror_ << "): " << H[1] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror_ << "): " << H[2] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror_ << "): " << H[3] << endl;
		
		matrix return_hessian(n_var, n_var);
		for (size_t row_index=0; row_index < n_var; ++row_index) {
			for (size_t col_index=0; col_index < n_var; ++col_index) {
				return_hessian.set(row_index, col_index, H[row_index*n_var + col_index]);
			}
		}
		return return_hessian;
	}
	
	
	int num_lower_ieq() {
		return inequality_constraints_lower_.size();
	}
	int num_upper_ieq() {
		return inequality_constraints_upper_.size();
	}
private:
protected:
	bool PRINT_;
	ASL *asl_;
	FILE *nl_;
	// real *J_;
	
	vector<double> X0_;
	fint *nerror_;
	vector<int> equality_constraints_;
	vector<int> inequality_constraints_lower_;
	vector<int> inequality_constraints_upper_;
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
		if (PRINT) jac_ci.print();
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
		for (size_t i=0; i<eq_items.size(); ++i) {
			eq_violation += abs(eq_items[i]);
		}
		double ieq_violation = 0.0;
		for (size_t i=0; i<ieq_items.size(); ++i) {
			// if (ieq_items[i] > 0){
			ieq_violation += bracket_plus(ieq_items[i]);
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
				num_variables_(nlp.num_primal() + 2*nlp.num_dual()), num_constraints_(nlp.num_dual()),
				num_nlp_variables_(nlp.num_primal()), num_nlp_constraints_eq_(nlp.num_dual_eq()), num_nlp_constraints_ieq_(nlp.num_dual_ieq()),
				hessian_(num_variables_, num_variables_), 
				jacobian_(num_constraints_, num_variables_), 
				gradient_(num_variables_),lower_bound_(num_variables_), upper_bound_(num_variables_), 
				jacobian_lower_bound_(num_constraints_), jacobian_upper_bound_(num_constraints_)
	 {
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		matrix hess = nlp_->lagrangian_hessian(iterate);
		matrix je = nlp_->constraints_equality_jacobian(iterate);
		matrix ji = nlp_->constraints_inequality_jacobian(iterate);
	
		int num_qp_variables = iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_);
		int num_qp_constraints = iterate.num_dual_eq_ + iterate.num_dual_ieq_;
		
		// cout << "num dual eq: " << iterate.num_dual_eq_ << endl;
		// cout << endl << "=======" << endl;
		// cout << je << endl;
		// cout << "=======" << endl;
		for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(eq_constraint_index,variables, je.get(eq_constraint_index,variables));
			}
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+eq_constraint_index, -1.0);
			jacobian_.set(eq_constraint_index,nlp_->num_primal()+iterate.num_dual_eq_+eq_constraint_index, +1.0);
		}
		
		// cout << "num dual eq: " << iterate.num_dual_ieq_ << endl;
		// cout << endl << "=======" << endl;
		// cout << ji << endl;
		// cout << "=======" << endl;	
		
		for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,variables, ji.get(ieq_constraint_index,variables));
				// jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,variables, ji.get(ieq_constraint_index,variables));
				// cout << "about to set"
			}
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+ieq_constraint_index, -1.0);
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index, nlp_->num_primal()+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+ieq_constraint_index, +1.0);
		}
		// cout << endl << "=======" << endl;
		// cout << jacobian_ << endl;
		// cout << "=======" << endl;
		
		for (size_t r=0; r<iterate.num_primal_; ++r) {
			for (size_t c=0; c<iterate.num_primal_; ++c) {
				hessian_.set(r,c, hess.get(r,c));
			}
		}
				
		// NLP gradient copied over
		for (size_t i=0; i<iterate.num_primal_; ++i)
			gradient_[i] = iterate.penalty_parameter_*grad_obj[i];
		// equality slacks both have positive:
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) {
			gradient_[iterate.num_primal_+i] = 1.0;
			gradient_[iterate.num_primal_+iterate.num_dual_eq_+i] = 1.0;
		}
		// inequality slacks only penalize the positive parts:
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) {
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+i] = 1.0;
			gradient_[iterate.num_primal_+2*iterate.num_dual_eq_+iterate.num_dual_ieq_+i] = 0.0;
		}
		// cout << endl << "gradient: " << gradient_ << endl;
		
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
		
		// print();
	}
	
	void inc_regularization(double reg) {
		for (size_t r=0; r<num_variables_; ++r) {
			for (size_t c=0; c<num_variables_; ++c) {
				double entry = hessian_.get(r,c);
				if (r == c)
					entry += reg;
				hessian_.set(r,c, entry);
			}
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
	int num_variables_, num_constraints_;
	int num_nlp_variables_, num_nlp_constraints_eq_, num_nlp_constraints_ieq_;
	matrix hessian_;
	matrix jacobian_;
	vector<double> gradient_;
	vector<double> lower_bound_;
	vector<double> upper_bound_;
	vector<double> jacobian_lower_bound_;
	vector<double> jacobian_upper_bound_;
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
		for (size_t stationarity_index=0; stationarity_index < stationarity.size(); ++stationarity_index) {
			// J(x)y
			for (size_t dual_eq_index=0; dual_eq_index < nlp_->num_dual_eq(); ++dual_eq_index){
				stationarity[rho_index] += jac_ce.get(dual_eq_index, rho_index)*eq_signflip*constraint_eq_dual_values[dual_eq_index];
			}
			
			// barJ(x)bary:
			for (size_t dual_ieq_index=0; dual_ieq_index < nlp_->num_dual_ieq(); ++dual_ieq_index){
				stationarity[rho_index] += jac_ci.get(dual_ieq_index, rho_index)*ieq_signflip*constraint_ieq_dual_values[dual_ieq_index];
			}
			
			rho[rho_index] = stationarity[stationarity_index];
			++rho_index;
		}
		// cout << "sta post: " << stationarity << endl << endl << endl;
		// cout << "sta: " << stationarity << endl << endl << endl;
		for (size_t constraint_eq_index=0; constraint_eq_index < nlp_->num_dual_eq(); ++constraint_eq_index) {
			rho[rho_index] = min(bracket_plus(constraint_eq_values[constraint_eq_index]), 1.0 - eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
			rho[rho_index] = min(bracket_minus(constraint_eq_values[constraint_eq_index]), 1.0 + eq_signflip*constraint_eq_dual_values[constraint_eq_index]);
			++rho_index;
		}
		
		for (size_t constraint_ieq_index=0; constraint_ieq_index < nlp_->num_dual_ieq(); ++constraint_ieq_index) {
			rho[rho_index] = min(bracket_plus(constraint_ieq_values[constraint_ieq_index]), 1.0 - ieq_signflip*constraint_ieq_dual_values[constraint_ieq_index]);
			++rho_index;
			// cout << "constraint_ieq_values[constraint_ieq_index=" << constraint_ieq_index << "]: " << constraint_ieq_values[constraint_ieq_index] << endl;
			// cout << "constraint_ieq_dual_values[constraint_ieq_index=" << constraint_ieq_index << "]: " << constraint_ieq_dual_values[constraint_ieq_index] << endl;
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
		bool PRINT=true;
		// cout << "OPERATOR FOR RESID FUNC: " << iterate.penalty_parameter_ << endl;
		
		vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f:
			stationarity[stationarity_index] *= iterate.penalty_parameter_;
		}
		
		vector<double> eq_con_values = nlp_->constraints_equality(iterate);
		vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
		
		return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, iterate.dual_eq_values_, iterate.dual_ieq_values_);
	}
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// first entry of rho(x,y,bary,mu):
		bool PRINT=true;
		vector<double> stationarity = nlp_->objective_gradient(iterate);
		
		for (size_t stationarity_index=0; stationarity_index < nlp_->num_primal(); ++stationarity_index) {
			// mu*grad f:
			stationarity[stationarity_index] *= iterate.penalty_parameter_;
			
			// H*d:
			for (size_t primal_index=0; primal_index < nlp_->num_primal(); ++primal_index) {
				stationarity[stationarity_index] += subproblem.hessian_.get(stationarity_index, primal_index)*step.primal_values_[primal_index];
			}
		}
		
		vector<double> eq_con_values = nlp_->constraints_equality(iterate);
		matrix jacobian_eq = nlp_->constraints_equality_jacobian(iterate);
		for (size_t eq_con_index=0; eq_con_index< nlp_->num_dual_eq(); ++eq_con_index) {
			// c_k + J_k^Td:
			for (size_t primal_index=0; primal_index < nlp_->num_primal(); ++primal_index) {
				eq_con_values[eq_con_index] += jacobian_eq.get(eq_con_index, primal_index)*step.primal_values_[primal_index];
			}
		}
		
		vector<double> ieq_con_values = nlp_->constraints_inequality(iterate);
		matrix jacobian_ieq = nlp_->constraints_inequality_jacobian(iterate);
		for (size_t ieq_con_index=0; ieq_con_index< nlp_->num_dual_ieq(); ++ieq_con_index) {
			// barc_k + barJ_k^Td:
			for (size_t primal_index=0; primal_index < nlp_->num_primal(); ++primal_index) {
				ieq_con_values[ieq_con_index] += jacobian_ieq.get(ieq_con_index, primal_index)*step.primal_values_[primal_index];
			}
		}
		
		return resid_helper(iterate, stationarity, eq_con_values, ieq_con_values, step.dual_eq_values_, step.dual_ieq_values_);
	}
protected:
};

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual()), first_(true)  {}
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem) {
		// qpOASES::SQProblem example_(nlp_->num_primal() + 2*nlp_->num_dual(), nlp_->num_dual());
		/* Setting up QProblem object. */
		example_.setPrintLevel(qpOASES::PL_NONE);
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
		opt.setToReliable();
		// opt.enableRegularisation = qpOASES::BooleanType(true);

		example_.setOptions(opt);
		cout.flush();
		cerr.flush();
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
		
		// real_t xOpt[6];
		vector<double> full_primal(nlp_->num_primal() + 2*nlp_->num_dual());
		example_.getPrimalSolution( &full_primal[0] );
		for (size_t primal_index=0; primal_index<nlp_->num_primal(); ++primal_index)
			step.primal_values_[primal_index] = full_primal[primal_index];
		// cout << "got primal" << endl;
		
		
		// cout << "primal: " << step.primal_values_ << endl;
		// 	- every QP variable has a multiplier:
		//		-- iterate.num_primal_: primal variables in NLP
		// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
		//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
		//	- every constraint has a multiplier:
		//		-- iterate.num_dual_eq_
		//		-- iterate.num_dual_ieq_
		// so in total: 
		vector<double> yOpt(subproblem.num_variables_ + subproblem.num_constraints_);
		assert(&yOpt[0] != NULL);
		// cout << "got temp dual" << endl;
		cerr.flush();
		cout.flush();
		example_.getDualSolution( &yOpt[0] );				
		
		step.status_ = ret;
		// to get eq con multipliers, read past NLP & slack variables:
		for (size_t i=0; i<subproblem.num_nlp_constraints_eq_; ++i) {
			step.dual_eq_values_[i] = -yOpt[subproblem.num_variables_+i];
		}
		// to get inequality constraint multipliers, read past NLP & slack variables and eq con multipliers:
		for (size_t i=0; i<subproblem.num_nlp_constraints_ieq_; ++i) {
			step.dual_ieq_values_[i] = -yOpt[subproblem.num_variables_+subproblem.num_nlp_constraints_eq_+i];
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
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			dot_product += gradient[i]*step.primal_values_[i];
		}
		return iterate.penalty_parameter_*(f + dot_product) + constraint_violation_func_(iterate, step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};

class LinearDecreaseFunction : public FunctionWithNLPState {
public:
	LinearDecreaseFunction(Nlp &nlp) : FunctionWithNLPState(nlp), constraint_violation_func_(nlp) {
		
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) const {
		bool PRINT=false;
		// double f = nlp_->objective(iterate);
		vector<double> gradient = nlp_->objective_gradient(iterate);
		// dotproduct(gradient, )
		double dot_product=0.0;
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			dot_product += gradient[i]*step.primal_values_[i];
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
		iSQOIterate it_and_step(2,1,1);
		bool PRINT = false;
		double alpha = 1.0;
		double penfunc_start = penalty_func_(iterate);
		double linear_reduction_start = linear_decrease_func_(iterate, step);
		
		double eta = 1e-8;
		
		for (int alpha_cuts = 0; alpha_cuts < 20; alpha_cuts++) {
			it_and_step.update(iterate, alpha, step);
			double penfunc_step = penalty_func_(it_and_step);
			if (PRINT) cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << endl;
			if (penfunc_step <= penfunc_start - eta*alpha*linear_reduction_start) {
				break;
			}
			alpha = 0.5*alpha;
		}
		return alpha;
	}
private:
protected:
	PenaltyFunction penalty_func_;
	LinearDecreaseFunction linear_decrease_func_;
};
class HessianShifter : public FunctionWithNLPState {
public:
	HessianShifter(Nlp &nlp) : FunctionWithNLPState(nlp), solve_qp_(nlp) {
		
	}
	iSQOStep operator()(const iSQOIterate &iterate, iSQOQuadraticSubproblem &subproblem) {
		
		iSQOStep return_step = solve_qp_(subproblem);
		double current_regularization = 1.0e1;
		int current_regularization_steps = 0;
		while ((return_step.status_ != 0) || (return_step.x_norm() > 1e9)) {
			current_regularization *= 10;
			// cerr << "--> retrying with " << current_regularization << endl;
			subproblem.inc_regularization(current_regularization);
			// cout << "subproblem: " <<endl;
			// subproblem.print();
			return_step = solve_qp_(subproblem);
			++current_regularization_steps;
			
			if (current_regularization_steps == 20){
				cout << "FAILURE!" << endl;
				break;
			}
		}
		if (current_regularization_steps == 0) current_regularization = 0.0;
		// if (current_regularization_steps == 20) 
		return return_step;
	}
private:
protected:
	SolveQuadraticProgram solve_qp_;
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
		printf("-------------------------------------------|");
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
		printf("-------------------------------------------|");
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
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, string step_type, double step_mix, double alpha) {
		// step, feas_iter, penalty_iter, alpha
		//  -> params: step, feas_iter, penalty_iter, alpha
		printf(output_format_post_,
				step_type.c_str(),
				step_mix, combination_step.x_norm(), constraint_violation_func_(penalty_iterate) - constraint_violation_func_(penalty_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
				alpha
					);
	}
protected:
	// Nlp *nlp_;
	ConstraintViolationFunction constraint_violation_func_;
	LinearDecreaseFunction linear_decrease_func_;
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
const char TextOutput::output_desc_post_[] = " TT    ||d||   FeasRed     PenRed |    alpha\n";
const char TextOutput::output_format_pre_[] = " %3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
const char TextOutput::output_format_subprob_[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |";
const char TextOutput::output_format_post_[] = " %s %9.2e %9.2e %9.2e %9.2e |%9.2e\n";


bool assert_close(double val1, double val2, double tol) {
	assert(abs(val1) <= abs(val2) + tol*abs(val1));
	assert(abs(val2) <= abs(val1) + tol*abs(val2));
	assert(((val1 > tol) && (val2 > tol)) || ((val1 < -tol) && (val2 < -tol)) || (val1 == 0 && val2 == 0));
	return (abs(val1) <= abs(val2) + tol*abs(val1)) && (abs(val2) <= abs(val1) + tol*abs(val2));
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
	
	string problem_file("/Users/traviscj/optimization/cute_nl_nopresolve/hs008.nl");
	if (argc>1) {
		problem_file = string(argv[1]);
	}
	AmplNlp problem(problem_file);//.c_str());
	// Hs014 problem;
	iSQOIterate penalty_iterate = problem.initial();
	penalty_iterate.penalty_parameter_ = 1.0e-1;

	iSQOIterate feasibility_iterate = problem.initial();
	feasibility_iterate.penalty_parameter_ = 0.0;
	
	// Utilities for NLP:
	PenaltyFunction penfunc(problem);
	ResidualFunction residual_func(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraint_violation(problem);
	LineSearchFunction linesearch(problem);
	LinearDecreaseFunction linear_decrease_func(problem);
	HessianShifter hessian_shifter_func(problem);
	// Other
	TextOutput text_output(problem);
	
	double epsilon = 1e-1;

	text_output.start();
	int iter=-1;
	// ALGORITHM A // Step 1
	for (iter = 0; iter < 200; iter ++ ) {
		iSQOStep combination_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
		
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		
		if ((residual_func(penalty_iterate) < 1e-6) && (constraint_violation(penalty_iterate) < 1e-6)) {
			cout << endl << "Termination 2a" << endl;
			break;
		}
		// if ((residual_func(feasibility) < 1e-6) && (constraint_violation(penalty_iterate) > 1e-6)) {
			// cout << endl << "Termination 2b" << endl;
			// break;
		// }

		iSQOQuadraticSubproblem penalty_subproblem(problem, penalty_iterate);
		
		iSQOStep penalty_step = hessian_shifter_func(penalty_iterate, penalty_subproblem);
		// cout << "penalty step: " << endl;
		// penalty_step.print();
		
		// TODO: move this into functor.
		
		// cout << "===================="<< endl;
		// cout << endl << "penalty res: " << residual_func(penalty_iterate,penalty_subproblem,penalty_step) << endl;
		// cout << "===================="<< endl;
		// 
		iSQOQuadraticSubproblem feasibility_subproblem(problem, feasibility_iterate);
		// feasibility_subproblem.print();
		// break;
		// feasibility_step.status_ = 99;
		// for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
// 			feasibility_step.primal_values_[primal_index] = 0.0;
// 		}
// 		
		double linear_reduction_penalty = linear_decrease_func(penalty_iterate,penalty_step);
		double violation = constraint_violation(penalty_iterate);
		double step_mix = -1.0;
		bool PRINT=false;
		string steptype = "4a";
		if (PRINT) cout << endl << "SCENARIO A CHECK: " << linear_reduction_penalty << " >= epsilon*" << violation << " = " << epsilon*violation<< ": ";
		if (linear_reduction_penalty >= epsilon*violation) {
			if (PRINT) cout << "PASSES!";
			
			step_mix = 1.0;
			text_output.subproblem(0.0, feasibility_iterate, feasibility_subproblem, combination_step);
			// TODO convex combination function:
			for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
				combination_step.primal_values_[primal_index] = penalty_step.primal_values_[primal_index];
			}
			for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
				combination_step.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
			}
			for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
				combination_step.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
			}
			// combination_step = penalty_step;
		} else {
			// cout << "Couldn't get far enough, checking alternate steps: " << endl;
			// cout << ""
			iSQOStep feasibility_step = solve_qp(feasibility_subproblem);
			// cout << "delta l_k(d_k', mu=mu_k): " << linear_decrease_func(penalty_iterate,penalty_step) << endl;
			// cout.flush();
			// cerr.flush();
			// cout << "delta l_k(d_k'', mu=0): " << linear_decrease_func(feasibility_iterate,feasibility_step) << endl;
			
			if (linear_reduction_penalty >= epsilon*linear_decrease_func(feasibility_iterate,feasibility_step)) {
				step_mix = 1e0;
				steptype = "5a";
			} else {
				
				// double 
				step_mix = 1.0;
				
				// cerr << "HERE BE DRAGONS" << endl;
				for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
					combination_step.primal_values_[primal_index] = step_mix*penalty_step.primal_values_[primal_index] + (1.0-step_mix)*feasibility_step.primal_values_[primal_index];
				}
				// for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
// 					combination_step.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
// 				}
// 				for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
// 					combination_step.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
// 				}
				cerr << "delta l_k(lincomb, 0): "<< linear_decrease_func(feasibility_iterate,combination_step) 
					 << " ?>=? delta l_k(feas, 0): " << epsilon*linear_decrease_func(feasibility_iterate,feasibility_step)
					 << ": " << (linear_decrease_func(feasibility_iterate,combination_step)  >=epsilon*linear_decrease_func(feasibility_iterate,feasibility_step) + 1e-12)
					 << endl;
				int num_comb_reductions = 0;
				while (! (linear_decrease_func(feasibility_iterate,combination_step) >= epsilon*linear_decrease_func(feasibility_iterate,feasibility_step) + 1e-12)) {
					cerr << "reducing: " << endl;
					step_mix = .5*step_mix;
					cerr << "step_mix: " << step_mix << endl;
					for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
						combination_step.primal_values_[primal_index] = step_mix*penalty_step.primal_values_[primal_index] + (1.0-step_mix)*feasibility_step.primal_values_[primal_index];
					}
					// for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
					// 	combination_step.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
					// }
					// for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
					// 	combination_step.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
					// }
					++num_comb_reductions;
					if (num_comb_reductions == 20) break;
				}
				cerr << "final step_mix: " << step_mix << endl;
				
				steptype = "5b";
				
				
				// break;
			}
			// cerr << "FAILS!";
			// break;
			// TODO convex combination function:
			for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
				combination_step.primal_values_[primal_index] = step_mix*penalty_step.primal_values_[primal_index] + (1.0-step_mix)*feasibility_step.primal_values_[primal_index];
			}
			for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
				combination_step.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
			}
			for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
				combination_step.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
			}
			text_output.subproblem(0.0, feasibility_iterate, feasibility_subproblem, feasibility_step);
		}
		
		
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		double alpha = linesearch(penalty_iterate, combination_step);
		
		double current_regularization = 0.0;
		text_output.subproblem(current_regularization, penalty_iterate, penalty_subproblem, penalty_step);
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, steptype, step_mix, alpha);
		
		//////////////////////////
		// ALGORITHM A // STEP 6
		//////////////////////////
		penalty_iterate.update(penalty_iterate, alpha, combination_step);
		for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
			penalty_iterate.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
		}
		for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
			penalty_iterate.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
		}
		feasibility_iterate.update(penalty_iterate, alpha, combination_step);
		for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
			feasibility_iterate.dual_eq_values_[dual_eq_index] = penalty_step.dual_eq_values_[dual_eq_index];
		}
		for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
			feasibility_iterate.dual_ieq_values_[dual_ieq_index] = penalty_step.dual_ieq_values_[dual_ieq_index];
		}
	}
	cout << endl << endl;
	penalty_iterate.print();
	return 0;
}


