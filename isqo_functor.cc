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

class iSQOStep {
public:
	iSQOStep(int number_primal, int number_dual_eq, int number_dual_ieq) : 
				num_primal_(number_primal), 
				num_dual_eq_(number_dual_eq),
				num_dual_ieq_(number_dual_ieq),
				primal_values_(number_primal),
				dual_eq_values_(num_dual_eq_),
				dual_ieq_values_(num_dual_ieq_),
				status_(-42)
		{
			// penalty_parameter_ = 1e-1;
		}
		iSQOStep(const iSQOStep& s) {
			cerr << "copy being made!" << endl;
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
		// cout << "-- initializing Hs" << endl;
	}
	double objective(const iSQOIterate &iterate) {
		// cout << "--- objective" << endl;
		return (iterate.primal_values_[0] - 2)*(iterate.primal_values_[0] - 2) + (iterate.primal_values_[1]-1)*(iterate.primal_values_[1]-1);
	}
	vector<double> constraints_equality(const iSQOIterate &iterate) {
		vector<double> retval(iterate.num_dual_eq_);
		// double *retval = new double[iterate.num_primal_];
		// cout << "iterate: " << 
			// iterate.primal_values_[0] << ", " << iterate.primal_values_[1] << endl;
		// cout << "--- constraints" << endl;
		retval[0] = iterate.primal_values_[0] - 2*iterate.primal_values_[1] + 1;
		// cout << " c_e0 = " << retval[0] << endl;
		return retval;
	}
	vector<double> constraints_inequality(const iSQOIterate &iterate) {
		// double *retval = new double[2];
		vector<double> retval(iterate.num_dual_ieq_);
		// cout << "--- constraints" << endl;
		retval[0] = .25*iterate.primal_values_[0]*iterate.primal_values_[0] + iterate.primal_values_[1]*iterate.primal_values_[1] - 1;
		// retval[0] = 0.0;
		// retval[1] = 90.01;
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
		// cout << "lagrangian hessian: ieq: " << iterate.dual_ieq_values_[0] << endl;
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
	AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
		if (PRINT_) cout << "Constructing an AmplNlp" << endl;
		ASL *asl;
		asl = ASL_alloc(ASL_read_pfgh);
		cout << "Filename: " << stub_str << endl;
		FILE *nl = jac0dim(stub_str, 57);
		if (PRINT_) cout << "ran jac0dim: " << nl << endl;
			
		X0 = (real *)Malloc(n_var*sizeof(real));
	
		int status = pfgh_read(nl, ASL_return_read_err | ASL_findgroups);
	
		num_primal_ = n_var;
		X0_.resize(num_primal());
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			X0_[primal_index] = X0[primal_index];
		}
		
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
				// this scheme won't work, because if ampl_constraint_index = 0 is a double-sided bound, we just get duplicate zeros
				//  - could use fortran-style indices.
				//  - could.... use multiple vectors?
			} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY) {
				if (PRINT_) cout << "lower bound only";
				if (PRINT_) cout << "(push -ineq)";
				inequality_constraints_lower_.push_back(ampl_constraint_index);
			} else {
				cout << "unsupported";
			}
			if (PRINT_) cout << endl;
		}
		
		// TODO add variable bounds here, using similar logic to the above.
		
		fint *nerror = (fint *) Malloc(1*sizeof(fint));
		real *x = (real *)Malloc(num_primal()*sizeof(real));
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
			x[primal_index] = X0_[primal_index];

		real *con = (real *)Malloc(n_con*sizeof(real));
		conval(x, con, nerror);
		for (size_t ampl_constraint_index =0 ; ampl_constraint_index < n_con; ++ampl_constraint_index) {
			if (PRINT_) cout << "i=" << ampl_constraint_index << ": l=" << LUrhs[2*ampl_constraint_index] << " <= c(x_k)= " << con[ampl_constraint_index] << " <= u=" << LUrhs[2*ampl_constraint_index+1] << endl;
		}
		
		num_dual_eq_ = equality_constraints_.size();
		num_dual_ieq_ = inequality_constraints_lower_.size() + inequality_constraints_upper_.size();
		asl_ = asl;
	}
	
	iSQOIterate initial() {
		iSQOIterate initial(num_primal(), num_dual_eq(), num_dual_ieq());
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			initial.primal_values_[primal_index] = X0_[primal_index];
			cout << "initial_" << primal_index << ": " << X0_[primal_index] <<endl;
		}
		return initial;
	}
// 	void nonefunc() {
// 		ASL *asl = asl_;
// 		real *x = (real *)Malloc(n_var*sizeof(real));
// 		x[0] = 2.0;
// 		x[1] = 2.0;
// 	
// 		fint *nerror = (fint *) Malloc(1*sizeof(fint));
// 		
// 	
// 		// real *con = (real *)Malloc(n_con*sizeof(real));
// // 		conval(x, con, nerror);
// 		// cout << "constraint value(nerror = " << *nerror << "): " << con[0] << endl;
// 		// cout << "constraint value(nerror = " << *nerror << "): " << con[1] << endl;
// 		// 	
// 		// cout << "LUrhs: " << LUrhs[0] << endl;
// 		// cout << "LUrhs: " << LUrhs[1] << endl;
// 		// cout << "LUrhs: " << LUrhs[2] << endl;
// 		// cout << "LUrhs: " << LUrhs[3] << endl;
// 	
// 		// vector<int> equality_constraints;		// equality_constraints[i] is the index in AMPL jacobian where equality constraint i lies.
// 		// vector<int> inequality_constraints;		// inequality_constraints[i] is the index in AMPL jacobian where inequality constraint i lies.
// 		// vector<int> inequality_;
// 		// equality_constraints_ inequality_constraints_
// 		
// 		// 	
// 		// named_vector_print("EQUALITY: ", equality_constraints_);
// 		// named_vector_print("INEQUALITY: ", inequality_constraints_);
// 		// // cout << "]" << endl;
// 
// 	
// 		cout << "J(nerror = " << *nerror << "): " << J_[0] << endl;
// 		cout << "J(nerror = " << *nerror << "): " << J_[1] << endl;
// 		cout << "J(nerror = " << *nerror << "): " << J_[2] << endl;
// 		cout << "J(nerror = " << *nerror << "): " << J_[3] << endl;
// 	
// 		if (false) {
// 			real *test = (real *)Malloc(n_var*sizeof(real));
// 			congrd(0, x, test, nerror);
// 			cout << "test[0] = " << test[0] << endl;
// 			cout << "test[1] = " << test[1] << endl;
// 			congrd(1, x, test, nerror);
// 			cout << "test[0] = " << test[0] << endl;
// 			cout << "test[1] = " << test[1] << endl;
// 		}
// 	}
// 	
	// zeroth order NLP quantities:
	double objective(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		*nerror = 0;
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		real obj = objval(0, x, nerror);
		if (PRINT_) cout << "objective value(nerror = " << *nerror << "): " << obj << endl;
		return obj;
	}
	
	vector<double> constraints_equality(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		*nerror = 0;
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		real *con = (real *)Malloc(n_con*sizeof(real));
		// vector<double> con(n_con);
		conval(x, &con[0], nerror);
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
		fint *nerror = (fint *)Malloc(sizeof(fint));
		*nerror = 0;
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		real *con = (real *)Malloc(n_con*sizeof(real));
		conval(x, con, nerror);
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
		return inequality_constraint_evaluation;
	}
	
	// first order NLP quantities:
	vector<double> objective_gradient(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		*nerror = 0;
		// real *gradient = (real *)Malloc(n_var*sizeof(real));
		vector<double> return_gradient(n_var);
		
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		objgrd(0, x, &return_gradient[0], nerror);
		if (PRINT_) cout << "objective gradient(nerror = " << *nerror << "): " << return_gradient[0] << endl;
		if (PRINT_) cout << "objective gradient(nerror = " << *nerror << "): " << return_gradient[1] << endl;
		
		return return_gradient;
	}
	matrix constraints_equality_jacobian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
			x[primal_index] = iterate.primal_values_[primal_index];
		
		// J_ = (real *)Malloc(n_var*n_con*sizeof(real));
		// jacval(x,J_,nerror);
		// 
		// for (size_t jacobian_entry=0; jacobian_entry<n_var*n_con; ++jacobian_entry) {
		// 	if (PRINT_) cout << "J[i=" << jacobian_entry << "]: " << J_[jacobian_entry] << endl;
		// }
		// 
		matrix equality_constraint_jacobian(equality_constraints_.size(), n_var);
		// for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index<equality_constraints_.size(); ++isqo_eq_constraint_index){
		// 	if (PRINT_) cout << " - " << isqo_eq_constraint_index << ordinal(isqo_eq_constraint_index) << " equality_constraint maps to AMPL's "<< equality_constraints_[isqo_eq_constraint_index] << ordinal(equality_constraints_[isqo_eq_constraint_index])<< " constraint";
		// 	if (PRINT_) cout << ": " << "[";
		// 	for (size_t var_index=0; var_index < n_var; ++var_index) {
		// 		if (PRINT_) if (var_index != 0) cout << ", ";
		// 		if (PRINT_) cout << J_[n_con*var_index + equality_constraints_[isqo_eq_constraint_index]];
		// 		// equality_constraint_jacobian.set(isqo_eq_constraint_index, var_index, J_[n_var*equality_constraints[isqo_eq_constraint_index] + var_index]);
		// 		equality_constraint_jacobian.set(isqo_eq_constraint_index, var_index, J_[n_con*var_index + equality_constraints_[isqo_eq_constraint_index]]);
		// 		// TODO get around this copy.
		// 	}
		// 	if (PRINT_) cout << "]";
		// 	if (PRINT_) cout << endl;
		// }
		// size_t isqo_eq_constraint_index=0;
		if (PRINT_) cout << "equal: " << equality_constraints_ << endl;
		for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
			// cout << ""
			real *G = (real *)Malloc(n_var*sizeof(real));
			congrd(equality_constraints_[isqo_eq_constraint_index], x, G, nerror);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				// cout << "J_{" << isqo_ineq_lower_constraint_index << ", " << primal_index << "} = " << G[primal_index] << endl;
				equality_constraint_jacobian.set(isqo_eq_constraint_index, primal_index, G[primal_index]);
			}
			// ++isqo_ineq_constraint_index;
		}
		if (PRINT_) cout << "equality jacobian: [" << endl;
		if (PRINT_) equality_constraint_jacobian.print();
		if (PRINT_) cout << "]" << endl;
		return equality_constraint_jacobian;
	}
	matrix constraints_inequality_jacobian(const iSQOIterate &iterate) {
		// PRINT_=false;
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		
		real *x = (real *)Malloc(iterate.num_primal_*sizeof(real));
		for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index){
			x[primal_index] = iterate.primal_values_[primal_index];
			// cout << "x[" << primal_index << "]: " << x[primal_index] << endl;
		}
		J_ = (real *)Malloc(n_con*n_var*sizeof(real));
		jacval(x,J_,nerror);
		
		matrix inequality_constraint_jacobian(num_dual_ieq(), n_var);
		size_t isqo_ineq_constraint_index=0;
		if (PRINT_) cout << "lower: " << inequality_constraints_lower_ << endl;
		if (PRINT_) cout << "upper: " << inequality_constraints_upper_ << endl;
		for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
			// cout << ""
			real *G = (real *)Malloc(n_var*sizeof(real));
			congrd(inequality_constraints_lower_[isqo_ineq_lower_constraint_index], x, G, nerror);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				// cout << "J_{" << isqo_ineq_lower_constraint_index << ", " << primal_index << "} = " << G[primal_index] << endl;
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, -G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
			real *G = (real *)Malloc(n_var*sizeof(real));
			congrd(inequality_constraints_upper_[isqo_ineq_upper_constraint_index], x, G, nerror);
			for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
				// cout << "J_{" << isqo_ineq_upper_constraint_index << ", " << primal_index << "} = " << G[primal_index] << endl;
				inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, G[primal_index]);
			}
			++isqo_ineq_constraint_index;
		}
		if (PRINT_) cout << "======" << endl;
		if (PRINT_) cout << inequality_constraint_jacobian << endl;
		if (PRINT_) cout << "======" << endl;
		
		
		if (PRINT_) cout << "inequality: [" << endl;
		if (PRINT_) inequality_constraint_jacobian.print();
		if (PRINT_) cout << "]" << endl;
		// PRINT_=false;
		return inequality_constraint_jacobian;
	}
	
	// second order NLP quantities:
	matrix lagrangian_hessian(const iSQOIterate &iterate) {
		ASL *asl = asl_;
		fint *nerror = (fint *)Malloc(sizeof(fint));
		*nerror = 0;
		
		real *H = (real *)Malloc(n_var*n_var*sizeof(real));
		real *OW = (real *)Malloc(sizeof(real));
		OW[0] = iterate.penalty_parameter_;
		real *Y = (real *)Malloc(n_con*sizeof(real));
		
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
		fullhes(H, n_var, 0, OW, Y);
		if (iterate.penalty_parameter_ == 0) {
			// cout << "Oops!" << endl;
			real *H_f_only = (real *)Malloc(n_var*n_var*sizeof(real));
			real *OW_f_only = (real *)Malloc(sizeof(real));
			OW_f_only[0] = 1.0;
			real *Y_f_only = (real *)Malloc(n_con*sizeof(real));
			for (size_t dual_index=0; dual_index < n_con; ++dual_index) Y[dual_index] = 0.0;
			fullhes(H_f_only, n_var, 0, OW_f_only, Y_f_only);
			
			for (size_t row_index=0; row_index < n_var; ++row_index)
				for (size_t col_index=0; col_index < n_var; ++col_index) 
					H[row_index*n_var + col_index] -= H_f_only[row_index*n_var + col_index];
		}
		
		if (PRINT_) cout << "H(nerror = " << *nerror << "): " << H[0] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror << "): " << H[1] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror << "): " << H[2] << endl;
		if (PRINT_) cout << "H(nerror = " << *nerror << "): " << H[3] << endl;
		
		matrix return_hessian(n_var, n_var);
		for (size_t row_index=0; row_index < n_var; ++row_index) {
			for (size_t col_index=0; col_index < n_var; ++col_index) {
				return_hessian.set(row_index, col_index, H[row_index*n_var + col_index]);
			}
		}
		return return_hessian;
	}
private:
protected:
	bool PRINT_;
	ASL *asl_;
	real *J_;
	
	vector<double> X0_;
	
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
		
		if (PRINT) cout << "constraintviolationfunction e 0: " << con_values_eq[0] << endl;
		// J(x)*d
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		if (PRINT) jac_ce.print();
		for (size_t row=0; row<iterate.num_dual_eq_; ++row) {
			for (size_t col=0; col<iterate.num_primal_; ++col) {
				con_values_eq[row] += jac_ce.get(row,col)*step.primal_values_[col];
			}
		}
		if (PRINT) cout << "constraintviolationfunction e 1: " << con_values_eq[0] << endl;
		
		// \bar{J}(x)*d
		if (PRINT) cout << "constraintviolationfunction i 0: " << con_values_ieq[0] << endl;
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		if (PRINT) jac_ci.print();
		for (size_t row=0; row<iterate.num_dual_ieq_; ++row) {
			for (size_t col=0; col<iterate.num_primal_; ++col) {
				con_values_ieq[row] += jac_ci.get(row,col)*step.primal_values_[col];
			}
		}
		if (PRINT) cout << "constraintviolationfunction i 1: " << con_values_ieq[0] << endl;
		
		return two_vectors(con_values_eq, con_values_ieq);
	}
	
	double two_vectors(const vector<double> &eq_items, const vector<double> &ieq_items) const {
		double eq_violation = 0.0;
		for (size_t i=0; i<eq_items.size(); ++i) {
			eq_violation += abs(eq_items[i]);
		}
		double ieq_violation = 0.0;
		for (size_t i=0; i<ieq_items.size(); ++i) {
			if (ieq_items[i] > 0){
				ieq_violation += abs(ieq_items[i]);
			}
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
			// TODO fix identity matrix insertion here.
			jacobian_.set(eq_constraint_index,2, -1.0);
			jacobian_.set(eq_constraint_index,3, +1.0);
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
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,nlp_->num_primal()+2*ieq_constraint_index+0, -1.0);
			jacobian_.set(iterate.num_dual_eq_+ieq_constraint_index,nlp_->num_primal()+2*ieq_constraint_index+1, +1.0);
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
			jacobian_lower_bound_[iterate.num_dual_eq_+ieq_index] = -1e10;
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
	
	void print() {
		cout << endl;
		cout << "slack hessian: " << endl;
		hessian_.print();
		cout << "slack jacobian: " << endl;
		jacobian_.print();
		
		cout << "slack gradient: " << endl;
		for(size_t i=0; i<num_variables_; ++i) cout << gradient_[i] << " ";
		cout << endl;
		cout << "slack lb: " << endl;
		for(size_t i=0; i<num_variables_; ++i) cout << lower_bound_[i] << " ";
		cout << endl;
		cout << "slack ub: " << endl;
		for(size_t i=0; i<num_variables_; ++i) cout << upper_bound_[i] << " ";
		cout << endl;
		cout << "slack lbA: " << endl;
		for(size_t i=0; i<num_constraints_; ++i) cout << jacobian_lower_bound_[i] << " ";
		cout << endl;
		cout << "slack ubA: " << endl;
		for(size_t i=0; i<num_constraints_; ++i) cout << jacobian_upper_bound_[i] << " ";
		cout << endl;
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
	double operator()(const iSQOIterate &iterate) const {
		bool PRINT=false;
		// first entry of rho(x,y,bary,mu):
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		vector<double> rho_1(iterate.num_primal_);
		vector<double> rho_2(iterate.num_dual_eq_), rho_3(iterate.num_dual_eq_);
		vector<double> rho_4(iterate.num_dual_ieq_), rho_5(iterate.num_dual_ieq_);
		if (PRINT) cout << endl;
		if (PRINT) cout << "==========" << endl;
		if (PRINT) cout << "Je:";
		if (PRINT) jac_ce.print();
		if (PRINT) cout << "; Ji:";
		if (PRINT) jac_ci.print();
		if (PRINT) cout << "; ye:" << iterate.dual_eq_values_[0] << ", yi:" << iterate.dual_ieq_values_[0] << endl;
		if (PRINT) cout << "==========" << endl;
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			rho_1[i] = 0.0;
			rho_1[i] += iterate.penalty_parameter_*grad_obj[i];
			
			if (PRINT) cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; " << endl;
			for (size_t dual_eq_index=0; dual_eq_index < nlp_->num_dual_eq(); ++dual_eq_index) {
				if (PRINT) cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; "
							<< jac_ci.get(dual_eq_index,i)*iterate.dual_ieq_values_[dual_eq_index] << "; " << endl;
				rho_1[i] += jac_ce.get(dual_eq_index,i)*iterate.dual_eq_values_[dual_eq_index];
			}
			for (size_t dual_ieq_index=0; dual_ieq_index < nlp_->num_dual_ieq(); ++dual_ieq_index) {
				if (PRINT) cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; "
							<< jac_ci.get(dual_ieq_index,i)*iterate.dual_ieq_values_[dual_ieq_index] << "; " << endl;
				rho_1[i] += jac_ci.get(dual_ieq_index,i)*iterate.dual_ieq_values_[dual_ieq_index];
			}
			
			if (PRINT) cout << "rho_1[i=" << i << "]: " << rho_1[i] << endl;
		}
		// second & third entry of rho(iter)
		double stupid=-1;
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) {

			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho3debug_b: " << iterate.dual_eq_values_[i] << endl;
			rho_2[i] = min( max( con_values_eq[i], 0.0), 1.0-iterate.dual_eq_values_[i]);
			rho_3[i] = min( max(-con_values_eq[i], 0.0), 1.0+iterate.dual_eq_values_[i]);
			if (PRINT) cout << "rho_2[i=" << i << "]: " << rho_2[i] << endl;
			if (PRINT) cout << "rho_3[i=" << i << "]: " << rho_3[i] << endl;
		}
		// fourth & fifth entry of rho(iter)
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) {
			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho5debug_b: " << iterate.dual_ieq_values_[i] << endl;
			rho_4[i] = min( max( con_values_ieq[i], 0.0), 1.0-iterate.dual_ieq_values_[i]);
			rho_5[i] = min( max(-con_values_ieq[i], 0.0),     iterate.dual_ieq_values_[i]);
			if (PRINT) cout << "rho_4[i=" << i << "]: " << rho_4[i] << endl;
			if (PRINT) cout << "rho_5[i=" << i << "]: " << rho_5[i] << endl;
		}
		
		double total=0.0;
		for (size_t i=0; i<iterate.num_primal_; ++i) { total += rho_1[i]*rho_1[i]; }
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) { total += rho_2[i]*rho_2[i]+rho_3[i]*rho_3[i]; }
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) { total += rho_4[i]*rho_4[i]+rho_5[i]*rho_5[i]; }
		if (total > 1e10) {
			cerr << "Bad shit's going down...";
			exit(5);
		}
		return sqrt(total);
	}
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// first entry of rho(x,y,bary,mu):
		bool PRINT=false;
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		vector<double> rho_1(iterate.num_primal_);
		vector<double> rho_2(iterate.num_dual_eq_), rho_3(iterate.num_dual_eq_);
		vector<double> rho_4(iterate.num_dual_ieq_), rho_5(iterate.num_dual_ieq_);
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			rho_1[i] = 0.0;
			rho_1[i] += iterate.penalty_parameter_*grad_obj[i];
			for (size_t variable_index=0; variable_index < iterate.num_primal_; ++variable_index) {
				rho_1[i] += subproblem.hessian_.get(i, variable_index)*step.primal_values_[variable_index];
			}
			for (size_t dual_eq_index=0; dual_eq_index<iterate.num_dual_eq_; dual_eq_index++) {
				rho_1[i] += jac_ce.get(dual_eq_index,i)*step.dual_eq_values_[dual_eq_index];
			}
			for (size_t dual_ieq_index=0; dual_ieq_index<iterate.num_dual_ieq_; dual_ieq_index++) {
				rho_1[i] += jac_ci.get(dual_ieq_index,i)*step.dual_ieq_values_[dual_ieq_index];
			}
			size_t j=0; 	// TODO removeme
			if (PRINT) cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; "
				<< jac_ce.get(i,j)*iterate.dual_eq_values_[j] << "; "
				<< jac_ci.get(i,j)*iterate.dual_ieq_values_[j] << "; " << endl;
			if (PRINT) cout << "rho_1[i=" << i << "]: " << rho_1[i] << endl;
		}
		// second & third entry of rho(iter)
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		for (size_t con_eq_index=0; con_eq_index<iterate.num_dual_eq_; ++con_eq_index) {
			for (size_t variable_index=0 ; variable_index < iterate.num_primal_; ++variable_index){
				con_values_eq[con_eq_index] += subproblem.jacobian_.get(con_eq_index,variable_index)*step.primal_values_[variable_index];
			}
		
			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho3debug_b: " << iterate.dual_eq_values_[i] << endl;			
			rho_2[con_eq_index] = min( max( con_values_eq[con_eq_index], 0.0), 1.0-step.dual_eq_values_[con_eq_index]);
			rho_3[con_eq_index] = min( max(-con_values_eq[con_eq_index], 0.0), 1.0+step.dual_eq_values_[con_eq_index]);
			if (PRINT) cout << "rho_2[i=" << con_eq_index << "]: " << rho_2[con_eq_index] << endl;
			if (PRINT) cout << "rho_3[i=" << con_eq_index << "]: " << rho_3[con_eq_index] << endl;
		}
		// fourth & fifth entry of rho(iter)
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		for (size_t con_ieq_index=0; con_ieq_index<iterate.num_dual_ieq_; ++con_ieq_index) {
			if (PRINT) cout << "rho4debug_a: " << max( con_values_ieq[con_ieq_index], 0.0) << endl;
			if (PRINT) cout << "rho4debug_b: " << 1.0-iterate.dual_ieq_values_[con_ieq_index] << endl;
			if (PRINT) subproblem.jacobian_.print();
			if (PRINT) cout << " - " << subproblem.jacobian_.get(0,0) << " * " << step.primal_values_[0] << endl;
			if (PRINT) cout << " - " << subproblem.jacobian_.get(0,1) << " * " << step.primal_values_[1] << endl;
			if (PRINT) cout << " - " << subproblem.jacobian_.get(1,0) << " * " << step.primal_values_[0] << endl;
			if (PRINT) cout << " - " << subproblem.jacobian_.get(1,1) << " * " << step.primal_values_[1] << endl;
			for (size_t variable_index=0 ; variable_index < iterate.num_primal_; ++variable_index){
				if (PRINT) cout << "-- Term: " << subproblem.jacobian_.get(iterate.num_dual_eq_+con_ieq_index,variable_index)*step.primal_values_[variable_index] << endl;
				con_values_ieq[con_ieq_index] += subproblem.jacobian_.get(iterate.num_dual_eq_+con_ieq_index,variable_index)*step.primal_values_[variable_index];
			}
			if (PRINT) cout << "rho4debug_a: " << max( con_values_ieq[con_ieq_index], 0.0) << endl;
			if (PRINT) cout << "rho4debug_b: " << 1.0-iterate.dual_ieq_values_[con_ieq_index] << endl;
		
			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho5debug_b: " << iterate.dual_ieq_values_[i] << endl;
			rho_4[con_ieq_index] = min( max( con_values_ieq[con_ieq_index], 0.0), 1.0-step.dual_ieq_values_[con_ieq_index]);
			rho_5[con_ieq_index] = min( max(-con_values_ieq[con_ieq_index], 0.0),     step.dual_ieq_values_[con_ieq_index]);
			if (PRINT) cout << "rho_4[i=" << con_ieq_index << "]: " << rho_4[con_ieq_index] << endl;
			if (PRINT) cout << "rho_5[i=" << con_ieq_index << "]: " << rho_5[con_ieq_index] << endl;
		}
		
		double total=0.0;
		for (size_t i=0; i<iterate.num_primal_; ++i) { total += rho_1[i]*rho_1[i]; }
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) { total += rho_2[i]*rho_2[i]+rho_3[i]*rho_3[i]; }
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) { total += rho_4[i]*rho_4[i]+rho_5[i]*rho_5[i]; }
		if (total > 1e10) {
			cerr << "Bad shit's going down...";
			exit(5);
		}
		return sqrt(total);
	}
protected:
};

class SolveQuadraticProgram : public FunctionWithNLPState {
public:
	SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(nlp.num_primal() + 2*nlp.num_dual(), nlp.num_dual()), first_(true)  {}
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem) {
		/* Setting up QProblem object. */
		example_.setPrintLevel(qpOASES::PL_LOW);
		int nWSR = 2000;
		qpOASES::returnValue ret;
		// cout << subproblem.hessian_ <<endl;
		if (first_) {
			ret = example_.init( &subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR );
			first_ = false;
		} else{
			ret = example_.hotstart(&subproblem.hessian_.data_[0],
								 &subproblem.gradient_[0],
								 &subproblem.jacobian_.data_[0],
								 &subproblem.lower_bound_[0],
								 &subproblem.upper_bound_[0],
								 &subproblem.jacobian_lower_bound_[0],
								 &subproblem.jacobian_upper_bound_[0],
								 nWSR);
		}
		cerr.flush();
		cout.flush();

		size_t return_status = example_.getStatus();
		// getGlobalMessageHandler()->listAllMessages();
		if( ret != qpOASES::SUCCESSFUL_RETURN )
		        printf( "%s\n", qpOASES::getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
		
		iSQOStep step(subproblem.num_nlp_variables_,subproblem.num_nlp_constraints_eq_,subproblem.num_nlp_constraints_ieq_);
		// real_t xOpt[6];
		example_.getPrimalSolution( &step.primal_values_[0] );

		// 	- every QP variable has a multiplier:
		//		-- iterate.num_primal_: primal variables in NLP
		// 		-- 2*iterate.num_dual_eq_: positive & negative slacks for equalities
		//		-- 2*iterate.num_dual_ieq_: positive & negative slacks for inequalities
		//	- every constraint has a multiplier:
		//		-- iterate.num_dual_eq_
		//		-- iterate.num_dual_ieq_
		// so in total: 
		vector<double> yOpt(subproblem.num_variables_ + subproblem.num_constraints_);
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
		double f = nlp_->objective(iterate);
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
		printf("----------------------------------|");
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
		printf("----------------------------------|");
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
				residual_func_(penalty_iterate), residual_func_(penalty_iterate)
				);
	}
	void subproblem(double shift, const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) const {
		// shift, step, linear_decrease_func, residual_func
		//  -> params: shift, step
		//	-> state: linear_decrease_func, residual_func_
		printf(output_format_subprob_,
				0.0, step.status_, step.x_norm(),
				linear_decrease_func_(iterate, step), residual_func_(iterate, subproblem, step)
					);
	}
	void post(const iSQOIterate &feasibility_iterate, const iSQOIterate &penalty_iterate, const iSQOStep &combination_step, double alpha) {
		// step, feas_iter, penalty_iter, alpha
		//  -> params: step, feas_iter, penalty_iter, alpha
		printf(output_format_post_,
				combination_step.x_norm(), constraint_violation_func_(penalty_iterate) - constraint_violation_func_(penalty_iterate,combination_step), linear_decrease_func_(penalty_iterate, combination_step),
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
const char TextOutput::output_format_post_[] = " 4a %9.2e %9.2e %9.2e |%9.2e\n";


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
		
}
int main() {
	// double x[2];
	
	// x[0]=2.78;
	// x[1]=4.66;
	// iSQOIterate penalty_iterate(2,1,1);
	AmplNlp problem("/Users/traviscj/optimization/cute_nl_nopresolve/hs011.nl");
	iSQOIterate penalty_iterate = problem.initial();
	// penalty_iterate.primal_values_[0] = 1.125;
	// penalty_iterate.primal_values_[1] = .125;
	penalty_iterate.penalty_parameter_ = 1.0e-1;
	vector<double> eq_values = problem.constraints_equality(penalty_iterate);
	for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
		if (eq_values[dual_eq_index] < 0) {
			penalty_iterate.dual_eq_values_[dual_eq_index] = -1.0;
		} else if (eq_values[dual_eq_index] > 0) {
			penalty_iterate.dual_eq_values_[dual_eq_index] = +1.0;
		} else {
			penalty_iterate.dual_eq_values_[dual_eq_index] = 0.0;
		}
	}
	vector<double> ieq_values = problem.constraints_inequality(penalty_iterate);
	for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
		if (ieq_values[dual_ieq_index] < 0) {
			penalty_iterate.dual_ieq_values_[dual_ieq_index] = 0.0;
		} else if (ieq_values[dual_ieq_index] > 0) {
			penalty_iterate.dual_ieq_values_[dual_ieq_index] =+1.0;
		} else {
			penalty_iterate.dual_ieq_values_[dual_ieq_index] = 0.0;
		}
	}
	// penalty_iterate.primal_values_[0] = 2.0;
	// penalty_iterate.primal_values_[1] = 2.0;
	// penalty_iterate.dual_eq_values_[0] = -1.0*.1;
	// penalty_iterate.dual_ieq_values_[0] = 0*.1;
	// penalty_iterate.dual_eq_values_[0] = -1.0;
	// penalty_iterate.dual_ieq_values_[0] = 1.0;

	
	iSQOIterate feasibility_iterate = problem.initial();
	feasibility_iterate.penalty_parameter_ = 0.0;
	// vector<double> eq_values = problem.constraints_equality(penalty_iterate);
	for (size_t dual_eq_index=0; dual_eq_index<problem.num_dual_eq(); ++dual_eq_index) {
		if (eq_values[dual_eq_index] < 0) {
			feasibility_iterate.dual_eq_values_[dual_eq_index] = -1.0;
		} else if (eq_values[dual_eq_index] > 0) {
			feasibility_iterate.dual_eq_values_[dual_eq_index] = +1.0;
		} else {
			feasibility_iterate.dual_eq_values_[dual_eq_index] = 0.0;
		}
	}
	// vector<double> ieq_values = problem.constraints_inequality(penalty_iterate);
	for (size_t dual_ieq_index=0; dual_ieq_index<problem.num_dual_ieq(); ++dual_ieq_index) {
		if (ieq_values[dual_ieq_index] < 0) {
			feasibility_iterate.dual_ieq_values_[dual_ieq_index] = 0.0;
		} else if (ieq_values[dual_ieq_index] > 0) {
			feasibility_iterate.dual_ieq_values_[dual_ieq_index] =+1.0;
		} else {
			feasibility_iterate.dual_ieq_values_[dual_ieq_index] = 0.0;
		}
	}
	// feasibility_iterate.primal_values_[0] = 2.0;
	// feasibility_iterate.primal_values_[1] = 2.0;
	// feasibility_iterate.dual_eq_values_[0] = -1*.1;
	// feasibility_iterate.dual_ieq_values_[0] = 0*.1;
	// feasibility_iterate.dual_eq_values_[0] = -1;
	// feasibility_iterate.dual_ieq_values_[0] = 0;
	
	// iSQOIterate final(2,1,1);
	// final.penalty_parameter_ = 1.0e-1;
	// final.primal_values_[0] = .822875656;
	// final.primal_values_[1] = .911437828;
	// final.dual_eq_values_[0] = 1.594493103554523e+00*.1;
	// final.dual_ieq_values_[0] = 1.846589027861980e+00*.1;
	// iSQOIterate final = penalty_iterate;
	
	// penalty_iterate = final;
	// penalty_iterate.print();
	
	// NLP object:
	// Hs014 problem;
	
	// Utilities for NLP:
	PenaltyFunction penfunc(problem);
	ResidualFunction residual_func(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraint_violation(problem);
	LineSearchFunction linesearch(problem);
	LinearDecreaseFunction linear_decrease_func(problem);
	// Other
	TextOutput text_output(problem);
	
	iSQOStep penalty_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
	iSQOStep feasibility_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());
	iSQOStep combination_step(problem.num_primal(),problem.num_dual_eq(),problem.num_dual_ieq());

	double epsilon = 1e-1;


	text_output.start();
	int iter=-1;
	// ALGORITHM A // Step 1
	for (iter = 0; iter < 35; iter ++ ) {
		// cout << "===================="<< endl;
		// cout << "residual: " << residual_func(penalty_iterate) << endl;
		// cout << "===================="<< endl;
		text_output.pre(iter, feasibility_iterate, penalty_iterate);
		if ((residual_func(penalty_iterate) < 1e-6) && (constraint_violation(penalty_iterate) < 1e-6)) {
			cout << endl << "Termination 2a" << endl;
			break;
		}

		iSQOQuadraticSubproblem penalty_subproblem(problem, penalty_iterate);
		penalty_step.status_ = 99;
		for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
			penalty_step.primal_values_[primal_index] = 0.0;
		}
		
		penalty_step = solve_qp(penalty_subproblem);
		
		// cout << "===================="<< endl;
		// cout << endl << "penalty res: " << residual_func(penalty_iterate,penalty_subproblem,penalty_step) << endl;
		// cout << "===================="<< endl;
		// 
		iSQOQuadraticSubproblem feasibility_subproblem(problem, feasibility_iterate);
		feasibility_step.status_ = 99;
		for (size_t primal_index=0; primal_index<problem.num_primal(); ++primal_index) {
			feasibility_step.primal_values_[primal_index] = 0.0;
		}
		
		double linear_reduction_penalty = linear_decrease_func(penalty_iterate,penalty_step);
		double violation = constraint_violation(penalty_iterate);
		double step_mix = -1.0;
		bool PRINT=false;
		if (PRINT) cout << endl << "SCENARIO A CHECK: " << linear_reduction_penalty << " >= epsilon*" << violation << " = " << epsilon*violation<< ": ";
		if (linear_reduction_penalty >= epsilon*violation) {
			if (PRINT) cout << "PASSES!";
			// combination_step = penalty_step;
			step_mix = 1.0;
		} else {
			cerr << "FAILS!";
			break;
		}
		
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
		
		//////////////////////////
		// ALGORITHM A // STEP 5
		//////////////////////////
		double alpha = linesearch(penalty_iterate, combination_step);
		
		text_output.subproblem(0.0, feasibility_iterate, feasibility_subproblem, feasibility_step);
		text_output.subproblem(0.0, penalty_iterate, penalty_subproblem, penalty_step);
		text_output.post(feasibility_iterate, penalty_iterate, combination_step, alpha);
		
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


