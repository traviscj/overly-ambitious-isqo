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
	double get(size_t r, size_t c) {
		return data_[columns_*r + c];
	}
	// double
	void print() {
		for (size_t r=0; r<rows_; r++) {
			for (size_t c=0; c<columns_; c++) {
				cout << " " << data_[columns_*r + c];
			}
			cout << endl;
		}
	}
	vector<double> multiply(const vector<double> x) {
		// for (size_t i=0; i<)
		vector<double> retval(3);
		retval[0] = 1;
		retval[1] = 2;
		retval[2] = 3;
		return retval;
	}
private:
protected:
	int rows_, columns_;
	vector<double> data_;
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

class NlpResidualFunction : public FunctionWithNLPState {
public:
	NlpResidualFunction(Nlp &nlp) : 
			FunctionWithNLPState(nlp) 
			// constraint_violation_func_(nlp),
			// penalty_parameter_(1.0) 
			{
		// cout << "-- Initializing nlp residual function." << endl;
	}
	double operator()(const iSQOIterate &iterate) {
		// first entry of rho(x,y,bary,mu):
		// vector<double> grad_obj(iterate.num_primal_);
		// vector<double> grad_ce(iterate.num_primal_);
		// vector<double> grad_ci(iterate.num_primal_);
		// matrix jac_ce(1,2);
		// matrix jac_ci(1,2);
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		// cout << "g:" << grad_obj[0] << ", " << grad_obj[1] << endl;
		matrix jac_ce = nlp_->constraints_equality_jacobian(iterate);
		// cout << "Je:" << jac_ce[0] << ", " << jac_ce[1] << endl;
		matrix jac_ci = nlp_->constraints_inequality_jacobian(iterate);
		// cout << "Ji:" << grad_ci[0] << ", " << grad_ci[1] << endl;
		double rho_1[iterate.num_primal_], rho_2[iterate.num_dual_eq_], rho_3[iterate.num_dual_eq_];
		double rho_4[iterate.num_dual_ieq_], rho_5[iterate.num_dual_ieq_];
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
protected:
};

USING_NAMESPACE_QPOASES
class GenerateQuadraticProgram : public FunctionWithNLPState {
public:
	GenerateQuadraticProgram(Nlp &nlp) : FunctionWithNLPState(nlp), example_(6,2), first_(true)  {}
	
	iSQOStep generateQP(const iSQOIterate &iterate) {
		// vector<double> grad_obj(iterate.num_primal_);
		vector<double> grad_obj = nlp_->objective_gradient(iterate);
		matrix hess = nlp_->lagrangian_hessian(iterate);
		// hess.print();
		matrix je = nlp_->constraints_equality_jacobian(iterate);
		// je.print();
		matrix ji = nlp_->constraints_inequality_jacobian(iterate);
		// ji.print();
	
		// TODO major cleanup needed here.
		matrix slack_jac(iterate.num_dual_eq_ + iterate.num_dual_ieq_, iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_));
		for (size_t eq_constraint_index=0; eq_constraint_index < iterate.num_dual_eq_; ++eq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				slack_jac.set(eq_constraint_index,variables, je.get(eq_constraint_index,variables));
			}
			slack_jac.set(eq_constraint_index,2, -1.0);
			slack_jac.set(eq_constraint_index,3, +1.0);
		}
		for (size_t ieq_constraint_index=0; ieq_constraint_index < iterate.num_dual_ieq_; ++ieq_constraint_index) {
			for (size_t variables=0; variables<iterate.num_primal_; ++variables) {
				slack_jac.set(iterate.num_dual_ieq_+ieq_constraint_index,variables, ji.get(ieq_constraint_index,variables));
			}
			// slack_jac.set(iterate.num_dual_ieq_+ieq_constraint_index,1, ji.get(ieq_constraint_index,1));
			slack_jac.set(iterate.num_dual_ieq_+ieq_constraint_index,2*iterate.num_dual_ieq_+2, -1.0);
			slack_jac.set(iterate.num_dual_ieq_+ieq_constraint_index,2*iterate.num_dual_ieq_+3, +1.0);
		}
		
		// cout << "slack jacobian:" << endl;
		// slack_jac.print();
	
		matrix slack_hess(iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_),iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_));
		for (size_t r=0; r<iterate.num_primal_; ++r) {
			for (size_t c=0; c<iterate.num_primal_; ++c) {
				slack_hess.set(r,c, hess.get(r,c));
			}
		}
		// cout << "Slack Hessian: " << endl;
		// slack_hess.print();
	
		vector<double> lower_bound_jac(iterate.num_dual_eq_ + iterate.num_dual_ieq_);
		vector<double> upper_bound_jac(iterate.num_dual_eq_ + iterate.num_dual_ieq_);
		vector<double> con_values_eq=nlp_->constraints_equality(iterate);
		vector<double> con_values_ieq=nlp_->constraints_inequality(iterate);
		
		for (size_t eq_index=0; eq_index<iterate.num_dual_eq_; ++eq_index) {
			lower_bound_jac[eq_index] = -con_values_eq[eq_index];
			upper_bound_jac[eq_index] = -con_values_eq[eq_index];
		}
		for (size_t ieq_index=0; ieq_index<iterate.num_dual_eq_; ++ieq_index) {
			lower_bound_jac[iterate.num_dual_eq_+ieq_index] = -1e10;
			upper_bound_jac[iterate.num_dual_eq_+ieq_index] = -con_values_ieq[ieq_index];
		}
	
	
		vector<double> lower_bound_var(iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_));
		vector<double> upper_bound_var(iterate.num_primal_ + 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_));
		for (size_t variable_index=0; variable_index < iterate.num_primal_; ++variable_index) {
			lower_bound_var[variable_index] = -1e10;
			upper_bound_var[variable_index] = +1e10;
		}
		for (size_t variable_index=0; variable_index < 2*(iterate.num_dual_eq_ + iterate.num_dual_ieq_); ++variable_index) {
			lower_bound_var[iterate.num_primal_+variable_index] = 0.0;
			upper_bound_var[iterate.num_primal_+variable_index] = 1e10;
			// cout << "lb[2+i="<<variable_index<< "]: " << lower_bound_var[iterate.num_primal_+variable_index] 
			// 	<< "ub[2+i="<<variable_index<< "]: " << upper_bound_var[iterate.num_primal_+variable_index] 
			// 	<< endl;
		}
		
		

		/* Setup data of first QP. */
		// for (size_t i=2; i<6; ++i) slack_hess.set(i,i,1e-8 + slack_hess.get(i,i));
		real_t H[6*6];// = { 1.0, 0.0, 0.0, 0.5 };
		for (size_t r=0; r<6; ++r) for(size_t c=0; c<6; ++c) H[6*r+c] = slack_hess.get(r,c);
		
		real_t A[2*6];// = { 1.0, 1.0 };
		for (size_t r=0; r<2; ++r) for(size_t c=0; c<6; ++c) A[6*r+c] = slack_jac.get(r,c);
		
		// for (int i=0; i<12; i++) cout << "A[i=" << i << "]: " << A[i] << endl;
		real_t g[6];
		g[0] = iterate.penalty_parameter_*grad_obj[0];
		g[1] = iterate.penalty_parameter_*grad_obj[1];
		g[2] = 1.0;
		g[3] = 1.0;
		g[4] = 1.0;
		g[5] = 0.0;
		
		
		real_t lb[6];
		for (size_t i=0; i<6; ++i) lb[i] = lower_bound_var[i];
		
		real_t ub[6];
		for (size_t i=0; i<6; ++i) ub[i] = upper_bound_var[i];
		
		
		real_t lbA[2];
		for (size_t i=0; i<2; ++i) lbA[i] = lower_bound_jac[i];
		
		real_t ubA[2];
		for (size_t i=0; i<2; ++i) ubA[i] = upper_bound_jac[i];
		
		bool PRINT = true;
		if (PRINT) {
			cout << endl;
			cout << "slack gradient: " << endl;
			for(size_t i=0; i<6; ++i) cout << g[i] << " ";
			cout << endl;
			cout << "slack lb: " << endl;
			for(size_t i=0; i<6; ++i) cout << lb[i] << " ";
			cout << endl;
			cout << "slack ub: " << endl;
			for(size_t i=0; i<6; ++i) cout << ub[i] << " ";
			cout << endl;
			cout << "slack lbA: " << endl;
			for(size_t i=0; i<2; ++i) cout << lbA[i] << " ";
			cout << endl;
			cout << "slack ubA: " << endl;
			for(size_t i=0; i<2; ++i) cout << ubA[i] << " ";
			cout << endl;
			
		}
		/* Setting up QProblem object. */
		// SQProblem example_( 6,2 );
		example_.setPrintLevel(PL_NONE);
		
		// cout << "Just about ready..." << endl;
		// cout << " - 1: " << (H != NULL) << endl;
		// cout << " - 2: " << (A != NULL) << endl;
		// cout << " - 3: " << (g != NULL) << endl;
		// cout << " - 4: " << (lb != NULL) << endl;
		// cout << " - 5: " << (ub != NULL) << endl;
		// cout << " - 6: " << (lbA != NULL) << endl;
		// cout << " - 7: " << (ubA != NULL) << endl;

		
		/* Solve first QP. */
		int nWSR = 2000;
		if (first_) {
			example_.init( H,g,A,lb,ub,lbA,ubA, nWSR );
			first_ = false;
		} else{
			example_.hotstart(H,g,A,lb,ub,lbA,ubA, nWSR );
		}
		cerr.flush();
		cout.flush();

		size_t return_status = example_.getStatus();
		// if (return_status != 0) {
			// exit(10);
		// }
		/* Get and print solution of second QP. */
		real_t xOpt[6];
		// vector<double> primal_return(6);
		example_.getPrimalSolution( xOpt );
		// for (size_t i=0; i<2; ++i) return_primal[i] = xOpt[i];
		// cout << "norm of d: " << sqrt(return_primal[0]*return_primal[0] + return_primal[1]*return_primal[1]) << endl;
		// cout << "got primal" << endl;
		real_t *yOpt = new real_t[200];
		example_.getDualSolution( yOpt );
		
		if (PRINT) {
			cout << "status: " << return_status << endl;
			printf( "xOpt = [ %e, %e, %e, %e, %e, %e ];  objVal = %e\n", xOpt[0],xOpt[1],xOpt[2],xOpt[3],xOpt[4],xOpt[5],example_.getObjVal() );
			printf( "yOpt = [ %e, %e, %e, %e, %e, %e, %e, %e]\n", yOpt[0], yOpt[1], yOpt[2], yOpt[3], yOpt[4], yOpt[5], yOpt[6], yOpt[7]);
		}
		
		iSQOStep step(iterate.num_primal_,iterate.num_dual_eq_,iterate.num_dual_ieq_);
		step.status_ = return_status;
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			step.primal_values_[i] = xOpt[i];
		}
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) {
			step.dual_eq_values_[i] = -yOpt[iterate.num_primal_+i];
		}
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) {
			step.dual_ieq_values_[i] = -yOpt[iterate.num_primal_+iterate.num_dual_eq_+i];
		}
		return step;
	}
private:
protected:
	SQProblem example_;
	bool first_;
};

class LineSearchFunction : public FunctionWithNLPState {
public:
	LineSearchFunction(Nlp &nlp) : FunctionWithNLPState(nlp), penfunc_(nlp) {
		// cout << "initializing a linesearcher!" << endl;
	}
	double operator()(const iSQOIterate &iterate, const iSQOStep &step) {
		iSQOIterate it_and_step(2,1,1);
		double alpha = 1.0;
		double penfunc_start = penfunc_(iterate);
		
		for (int alpha_cuts = 0; alpha_cuts < 20; alpha_cuts++) {
			it_and_step.update(iterate, alpha, step);
			double penfunc_step = penfunc_(it_and_step);
			// cout << "alpha_cut: " << alpha_cuts << ": original: " << penfunc_start << "; new: " << penfunc_step << "; reduction: " << (penfunc_start-penfunc_step) << endl;
			if (penfunc_step-penfunc_start <= -1e-6) {
				break;
			}
			alpha = 0.5*alpha;
		}
		return alpha;
	}
private:
protected:
	PenaltyFunction penfunc_;
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
		
		return -iterate.penalty_parameter_*dot_product + constraint_violation_func_(iterate) - constraint_violation_func_(iterate,step);
	}
private:
protected:
	ConstraintViolationFunction constraint_violation_func_;
};
int main() {
	Hs014 problem;
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
	iterate.dual_eq_values_[0] = 0;
	iterate.dual_ieq_values_[0] = 1;
	
	iSQOIterate final(2,1,1);
	final.penalty_parameter_ = 1.0e-1;
	final.primal_values_[0] = .822875656;
	final.primal_values_[1] = .911437828;
	final.dual_eq_values_[0] = 1.846589027861980e+00;
	final.dual_ieq_values_[0] = 1.594493103554523e+00;
	
	     
    
	PenaltyFunction penfunc(problem);
	
	cout << "current penalty is: " << penfunc(sillystart) << endl;
	cout << "current penalty is: " << penfunc(iterate) << endl;
	cout << "current penalty is: " << penfunc(final) << endl;
	
	NlpResidualFunction residual_func(problem);
	cout << "residual @ final: " << residual_func(final) << endl;
	
	// Utilities for NLP:
	GenerateQuadraticProgram qp(problem);
	ConstraintViolationFunction constraintviolation(problem);
	LineSearchFunction linesearch(problem);
	
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

		step = qp.generateQP(iterate);		
		double alpha = linesearch(iterate, step);
		iterate.update(iterate, alpha, step);
		
		LinearDecreaseFunction linear_decrease_func(problem);
		printf(output_format_post,
				-42.1, -42, step.x_norm(),
				linear_decrease_func(iterate, step), -42.1,
				alpha
					);
		if (residual_func(iterate) < 1e-8) {
			break;
		}
	}
	printf(output_format_pre, 
			iter, 
			problem.objective(iterate), constraintviolation(iterate),
			iterate.penalty_parameter_, penfunc(iterate),
			-residual_func(iterate), residual_func(iterate)
			);
	printf("\n");
	
	return 0;
}