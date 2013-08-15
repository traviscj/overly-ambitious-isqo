#include <iostream>
#include <vector>
// #include <boost/array.hpp>
// #include <array>
#include <cmath>
#include <qpOASES.hpp>

using namespace std;

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
	double x_norm() {
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
	double x_norm() {
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
		for (size_t dual_eq_index=0; dual_eq_index < iterate.num_dual_eq_; ++dual_eq_index) {
			dual_eq_values_[dual_eq_index] = step.dual_eq_values_[dual_eq_index];
		}
		for (size_t dual_ieq_index=0; dual_ieq_index < iterate.num_dual_ieq_; ++dual_ieq_index) {
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
				cout << " " << data_[columns_*r + c];
			}
			cout << endl;
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
class Nlp {
	// this class implements 
public:
	Nlp() {
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
protected:
};

// TODO: implement this using AMPL instead.
// TODO: Then, implement it using sparse AMPL instead...
class Hs014 : public Nlp {
public:
	Hs014() {
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
		hess.set(0,0,iterate.penalty_parameter_*2.0 + 0.5*iterate.dual_ieq_values_[0]/iterate.penalty_parameter_);
		hess.set(1,1,iterate.penalty_parameter_*2.0 + 2.0*iterate.dual_ieq_values_[1]/iterate.penalty_parameter_);
		return hess;
	}
protected:
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
	double operator()(const iSQOIterate &iterate) {
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		
		return two_vectors(con_values_eq, con_values_ieq);
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
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
	
	double two_vectors(const vector<double> &eq_items, const vector<double> &ieq_items) {
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
	double operator()(const iSQOIterate &iterate) {
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
				num_variables_(6), num_constraints_(2),
				num_nlp_variables_(2), num_nlp_constraints_eq_(1), num_nlp_constraints_ieq_(1),
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
		
	
		for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(eq_constraint_index,variables, je.get(eq_constraint_index,variables));
			}
			jacobian_.set(eq_constraint_index,2, -1.0);
			jacobian_.set(eq_constraint_index,3, +1.0);
		}
		for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				jacobian_.set(iterate.num_dual_ieq_+ieq_constraint_index,variables, ji.get(ieq_constraint_index,variables));
			}
			jacobian_.set(iterate.num_dual_ieq_+ieq_constraint_index,2*iterate.num_dual_ieq_+2, -1.0);
			jacobian_.set(iterate.num_dual_ieq_+ieq_constraint_index,2*iterate.num_dual_ieq_+3, +1.0);
		}
		
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
			
		vector<double> con_values_eq=nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);
		
		for (size_t eq_index=0; eq_index<iterate.num_dual_eq_; ++eq_index) {
			jacobian_lower_bound_[eq_index] = -con_values_eq[eq_index];
			jacobian_upper_bound_[eq_index] = -con_values_eq[eq_index];
		}
		for (size_t ieq_index=0; ieq_index<iterate.num_dual_eq_; ++ieq_index) {
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
	ResidualFunction(Nlp &nlp) : 
			FunctionWithNLPState(nlp) 
			// constraint_violation_func_(nlp),
			// penalty_parameter_(1.0) 
			{
		// cout << "-- Initializing nlp residual function." << endl;
	}
	double operator()(const iSQOIterate &iterate) {
		// first entry of rho(x,y,bary,mu):
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		vector<double> rho_1(iterate.num_primal_);
		vector<double> rho_2(iterate.num_dual_eq_), rho_3(iterate.num_dual_eq_);
		vector<double> rho_4(iterate.num_dual_ieq_), rho_5(iterate.num_dual_ieq_);
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			rho_1[i] = 0.0;
			rho_1[i] += iterate.penalty_parameter_*grad_obj[i];
			size_t j=0; // TODO: make this a for loop instead...
			
			rho_1[i] += jac_ce.get(j,i)*iterate.dual_eq_values_[j];
			rho_1[i] += jac_ci.get(j,i)*iterate.dual_ieq_values_[j];
			// cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; "
				// << jac_ce.get(i,j)*iterate.dual_eq_values_[j] << "; "
				// << jac_ci.get(i,j)*iterate.dual_ieq_values_[j] << "; " << endl;
			// cout << "rho_1[i=" << i << "]: " << rho_1[i] << endl;
		}
		// second & third entry of rho(iter)
		vector<double> con_values_eq = nlp_->constraints_equality(iterate);
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) {

			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho3debug_b: " << iterate.dual_eq_values_[i] << endl;
			rho_2[i] = min( max( con_values_eq[i], 0.0), 1.0-iterate.dual_eq_values_[i]);
			rho_3[i] = min( max(-con_values_eq[i], 0.0), 1.0+iterate.dual_eq_values_[i]);
			// cout << "rho_2[i=" << i << "]: " << rho_2[i] << endl;
			// cout << "rho_3[i=" << i << "]: " << rho_3[i] << endl;
		}
		// fourth & fifth entry of rho(iter)
		vector<double> con_values_ieq = nlp_->constraints_inequality(iterate);
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) {
			// cout << "rho5debug_a: " << max(-con_values_ieq[i], 0.0) << endl;
			// cout << "rho5debug_b: " << iterate.dual_ieq_values_[i] << endl;
			rho_4[i] = min( max( con_values_ieq[i], 0.0), 1.0-iterate.dual_ieq_values_[i]);
			rho_5[i] = min( max(-con_values_ieq[i], 0.0),     iterate.dual_ieq_values_[i]);
			// cout << "rho_4[i=" << i << "]: " << rho_4[i] << endl;
			// cout << "rho_5[i=" << i << "]: " << rho_5[i] << endl;
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
	double operator()(const iSQOIterate &iterate, const iSQOQuadraticSubproblem &subproblem, const iSQOStep &step) {
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
			// cout << "debug rho1_a: " << iterate.penalty_parameter_*grad_obj[i] << "; "
				// << jac_ce.get(i,j)*iterate.dual_eq_values_[j] << "; "
				// << jac_ci.get(i,j)*iterate.dual_ieq_values_[j] << "; " << endl;
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
	SolveQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(6,2), first_(true)  {}
	
	iSQOStep operator()(const iSQOQuadraticSubproblem &subproblem) {
		/* Setting up QProblem object. */
		example_.setPrintLevel(qpOASES::PL_LOW);
		int nWSR = 2000;
		qpOASES::returnValue ret;
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
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
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



int main() {
	// double x[2];
	
	// x[0]=2.78;
	// x[1]=4.66;
	iSQOIterate sillystart(2,1,1);
	sillystart.penalty_parameter_ = 1.0e-1;
	sillystart.primal_values_[0] = 2.78;
	sillystart.primal_values_[1] = 4.66;
	sillystart.dual_eq_values_[0] = 137.1354;
	sillystart.dual_ieq_values_[0] = 484.418;
	
	iSQOIterate iterate(2,1,1);
	iterate.penalty_parameter_ = 1.0e-1;
	iterate.primal_values_[0] = 2.0;
	iterate.primal_values_[1] = 2.0;
	iterate.dual_eq_values_[0] = -1;
	iterate.dual_ieq_values_[0] = 0;
	
	iSQOIterate final(2,1,1);
	final.penalty_parameter_ = 1.0e-1;
	final.primal_values_[0] = .822875656;
	final.primal_values_[1] = .911437828;
	final.dual_eq_values_[0] = 1.846589027861980e+00;
	final.dual_ieq_values_[0] = 1.594493103554523e+00;
	
	// NLP object:
	Hs014 problem;
	// Utilities for NLP:
	PenaltyFunction penfunc(problem);
	ResidualFunction residual_func(problem);
	SolveQuadraticProgram solve_qp(problem);
	ConstraintViolationFunction constraintviolation(problem);
	LineSearchFunction linesearch(problem);
	LinearDecreaseFunction linear_decrease_func(problem);
	
	// cout << "current penalty is: " << penfunc(sillystart) << endl;
	// cout << "current penalty is: " << penfunc(iterate) << endl;
	// cout << "current penalty is: " << penfunc(final) << endl;
	// cout << "residual @ final: " << residual_func(final) << endl;
	
	iSQOStep step(2,1,1);
	char output_desc_pre[] = " it  |      obj     infeas |       pen      merit |   feaskkt     penkkt &";
	char output_desc_post[] = "     shift  msg     ||d||     penred        res  | isqo | alpha\n";
	char output_format_pre[] = "%3d | %9.2e  %9.2e | %9.2e  %+9.2e | %9.2e  %9.2e &";
	char output_format_post[] = " %9.2e  %3d  %9.2e  %9.2e  %9.2e |      | %9.2e\n";
	
	printf("%s",output_desc_pre);
	printf("%s",output_desc_post);
	int iter=-1;
	for (iter = 0; iter < 15; iter ++ ) {
		printf(output_format_pre, 
				iter, 
				problem.objective(iterate), constraintviolation(iterate),
				iterate.penalty_parameter_, penfunc(iterate),
				-residual_func(iterate), residual_func(iterate)
				);
		if (residual_func(iterate) < 1e-6) {
			break;
		}

		// step = qp.generateQP(iterate);
		iSQOQuadraticSubproblem penalty_subproblem(problem, iterate);
		step = solve_qp(penalty_subproblem);
		
		// test scenarios...
		double epsilon = 1e-1;
		
		// TEST SCENARIO A
		/////////////////////////////////
		double linear_reduction = linear_decrease_func(iterate,step);
		double violation = constraintviolation(iterate);
		bool PRINT=false;
		if (PRINT) cout << endl << "SCENARIO A CHECK: " << linear_reduction << " >= epsilon*" << violation << " = " << epsilon*violation<< ": ";
		if (linear_reduction >= epsilon*violation) {
			if (PRINT) cout << "PASSES!";
		} else {
			if (PRINT) cout << "FAILS!";
		}
		if (PRINT) cout << endl;
		
		// TEST SCENARIO B
		/////////////////////////////////
		
		// TEST SCENARIO C
		/////////////////////////////////
		
		double alpha = linesearch(iterate, step);
		
		// cout << endl << "Testing residual: " <<  << endl;
		
		printf(output_format_post,
				0.0, step.status_, step.x_norm(),
				linear_decrease_func(iterate, step),residual_func(iterate, penalty_subproblem, step),
				alpha
					);
		
		// NOW update the iterate:
		iterate.update(iterate, alpha, step);
		
	}
	cout << endl << endl;
	
	return 0;
}