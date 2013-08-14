#include <iostream>
#include <cmath>
using namespace std;

class iSQOIterate {
public:
	iSQOIterate(int number_primal, int number_dual_eq, int number_dual_ieq) : 
			num_primal_(number_primal), 
			num_dual_eq_(number_dual_eq),
			num_dual_ieq_(number_dual_ieq)
	{
		primal_values_ = new double[num_primal_];
		dual_eq_values_ = new double[num_dual_eq_];
		dual_ieq_values_ = new double[num_dual_ieq_];
		penalty_parameter_ = 1e-1;
	}
// private:
// protected:
	int num_primal_;
	int num_dual_eq_;
	int num_dual_ieq_;
	
	double *primal_values_;
	double *dual_eq_values_;
	double *dual_ieq_values_;
	double penalty_parameter_;
};

class Nlp {
	// this class implements 
public:
	Nlp() { cout << "-- initializing Nlp" << endl; }
	
	// zeroth order NLP quantities:
	virtual double objective(const iSQOIterate &iterate) = 0;
	virtual double *constraints_equality(const iSQOIterate &iterate, double *retval) = 0;
	virtual double *constraints_inequality(const iSQOIterate &iterate, double *retval) = 0;
	
	// first order NLP quantities:
	virtual void objective_gradient(const iSQOIterate &iterate, double *ret_gradient) = 0;
	virtual void constraints_equality_jacobian(const iSQOIterate &iterate, double *ret_jacobian) = 0;
	virtual void constraints_inequality_jacobian(const iSQOIterate &iterate, double *ret_jacobian) = 0;
protected:
};

class Hs014 : public Nlp {
public:
	Hs014() { cout << "-- initializing Hs" << endl;}
	double objective(const iSQOIterate &iterate) {
		cout << "--- objective" << endl;
		return (iterate.primal_values_[0] - 2)*(iterate.primal_values_[0]-2) + (iterate.primal_values_[1]-1)*(iterate.primal_values_[1]-1);
	}
	double *constraints_inequality(const iSQOIterate &iterate, double *retval) {
		// double *retval = new double[2];
		cout << "--- constraints" << endl;
		retval[0] = .25*iterate.primal_values_[0]*iterate.primal_values_[0] + iterate.primal_values_[1]*iterate.primal_values_[1] - 1;
		// retval[1] = 90.01;
		return retval;
	}
	double *constraints_equality(const iSQOIterate &iterate, double *retval) {
		// double *retval = new double[2];
		cout << "iterate: " << 
			iterate.primal_values_[0] << ", " << iterate.primal_values_[1] << endl;
		cout << "--- constraints" << endl;
		retval[0] = iterate.primal_values_[0] - 2*iterate.primal_values_[1] + 1;
		cout << " c_e0 = " << retval[0] << endl;
		return retval;
	}
	
	void objective_gradient(const iSQOIterate &iterate, double *ret_gradient) {
		ret_gradient[0] = 2*(iterate.primal_values_[0] - 2);
		ret_gradient[1] = 2*(iterate.primal_values_[1] - 1);
	}
	void constraints_equality_jacobian(const iSQOIterate &iterate, double *ret_jacobian){
		ret_jacobian[0] = 2.0*.25*iterate.primal_values_[0];
		ret_jacobian[1] = 2.0*iterate.primal_values_[1];
	}
	void constraints_inequality_jacobian(const iSQOIterate &iterate, double *ret_jacobian){
		ret_jacobian[0] = 1.0;
		ret_jacobian[1] = -2.0;
	}
protected:
};

class FunctionWithNLPState {
public:
	FunctionWithNLPState(Nlp *nlp) {
		nlp_ = nlp;
	}
protected:
	Nlp *nlp_;
};

class ConstraintViolationFunction : public FunctionWithNLPState{
public:
	ConstraintViolationFunction(Nlp *nlp) : FunctionWithNLPState(nlp){
		cout << "-- Initializing l1 violation function." << endl;
	}
	double operator()(const iSQOIterate &iterate) {
		double con_values_eq[iterate.num_dual_eq_];
		double con_values_ieq[iterate.num_dual_ieq_];
		nlp_->constraints_equality(iterate, con_values_eq);
		nlp_->constraints_inequality(iterate, con_values_ieq);
		for (size_t i=0; i<iterate.num_dual_ieq_; i++) {
			if (con_values_ieq[i] <= 0){
				con_values_ieq[i] = 0;
			}
		}
		cout << "--- constraints: " << endl;
		cout << "---- eq: ";
		double eq_violation = 0.0;
		for (size_t i=0; i<iterate.num_dual_eq_; i++) {
			eq_violation += con_values_eq[i];
			cout << con_values_eq[i] << ", ";
		}
		cout << endl;
		cout << "---- ieq: ";
		double ieq_violation = 0.0;
		for (size_t i=0; i<iterate.num_dual_eq_; i++) {
			ieq_violation += con_values_ieq[i];
			cout << con_values_ieq[i] << endl;
		}
		return eq_violation + ieq_violation;
	}
protected:
};

class PenaltyFunction : public FunctionWithNLPState{
public:
	PenaltyFunction(Nlp *nlp) : 
			FunctionWithNLPState(nlp), 
			constraint_violation_func_(nlp),
			penalty_parameter_(1.0) {
		cout << "-- Initializing penalty function." << endl;
	}
	double operator()(const iSQOIterate &iterate) {
		cout << "--- Calling the PenaltyFunction Functor..." << endl;
		double f = nlp_->objective(iterate);
		cout << "--- objective: " << f << endl;
		return penalty_parameter_*f + constraint_violation_func_(iterate);
	}
protected:
	ConstraintViolationFunction constraint_violation_func_;
	double penalty_parameter_;
};

class NlpResidualFunction : public FunctionWithNLPState {
public:
	NlpResidualFunction(Nlp *nlp) : 
			FunctionWithNLPState(nlp) 
			// constraint_violation_func_(nlp),
			// penalty_parameter_(1.0) 
			{
		cout << "-- Initializing nlp residual function." << endl;
	}
	double operator()(const iSQOIterate &iterate) {
		// first entry of rho(x,y,bary,mu):
		double grad_obj[iterate.num_primal_], grad_ce[iterate.num_primal_], grad_ci[iterate.num_primal_];
		nlp_->objective_gradient(iterate, grad_obj);
		cout << "g:" << grad_obj[0] << ", " << grad_obj[1] << endl;
		nlp_->constraints_equality_jacobian(iterate, grad_ce);
		cout << "Je:" << grad_ce[0] << ", " << grad_ce[1] << endl;
		nlp_->constraints_inequality_jacobian(iterate, grad_ci);
		cout << "Ji:" << grad_ci[0] << ", " << grad_ci[1] << endl;
		double rho_1[iterate.num_primal_], rho_2[iterate.num_dual_eq_], rho_3[iterate.num_dual_eq_];
		double rho_4[iterate.num_dual_ieq_], rho_5[iterate.num_dual_ieq_];
		for (size_t i=0; i<iterate.num_primal_; ++i) {
			rho_1[i] = 0.0;
			rho_1[i] += iterate.penalty_parameter_*grad_obj[i];
			rho_1[i] += grad_ce[i]*iterate.dual_eq_values_[0];
			rho_1[i] += grad_ci[i]*iterate.dual_ieq_values_[0];
			cout << "rho_1[i=" << i << "]: " << rho_1[i] << endl;
		}
		// second & third entry of rho(iter)
		double con_values_eq[iterate.num_dual_eq_];
		nlp_->constraints_equality(iterate, con_values_eq);
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) {
			rho_2[i] = min( max( con_values_eq[i], 0.0), 1.0-iterate.dual_eq_values_[i]);
			rho_3[i] = min( max(-con_values_eq[i], 0.0), 1.0+iterate.dual_eq_values_[i]);
			cout << "rho_2[i=" << i << "]: " << rho_2[i] << endl;
			cout << "rho_3[i=" << i << "]: " << rho_3[i] << endl;
		}
		// fourth & fifth entry of rho(iter)
		double con_values_ieq[iterate.num_dual_ieq_];
		nlp_->constraints_inequality(iterate, con_values_ieq);
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) {
			rho_4[i] = min( max( con_values_ieq[i], 0.0), 1.0-iterate.dual_ieq_values_[i]);
			rho_5[i] = min( max(-con_values_ieq[i], 0.0),     iterate.dual_ieq_values_[i]);
			cout << "rho_4[i=" << i << "]: " << rho_4[i] << endl;
			cout << "rho_5[i=" << i << "]: " << rho_5[i] << endl;
		}
		
		double total=0.0;
		for (size_t i=0; i<iterate.num_primal_; ++i) { total += rho_1[i]*rho_1[i]; }
		for (size_t i=0; i<iterate.num_dual_eq_; ++i) { total += rho_2[i]*rho_2[i]+rho_3[i]*rho_3[i]; }
		for (size_t i=0; i<iterate.num_dual_ieq_; ++i) { total += rho_4[i]*rho_4[i]+rho_5[i]*rho_5[i]; }
		return sqrt(total);
	}
protected:
};

// class NlpUtilityFunctions {
// public:
// 	NlpUtilityFunctions(Nlp *nlp) {
// 		nlp_ = nlp;
// 	}
// 	
// protected:
// 	Nlp* nlp_;
// };

int main() {
	Hs014 problem;
	// double x[2];
	
	// x[0]=2.78;
	// x[1]=4.66;
	iSQOIterate sillystart(2,1,1);
	sillystart.primal_values_[0] = 2.78;
	sillystart.primal_values_[1] = 4.66;
	sillystart.dual_eq_values_[0] = 2.78;
	sillystart.dual_ieq_values_[0] = 4.66;
	
	iSQOIterate iterate(2,1,1);
	iterate.primal_values_[0] = 2.0;
	iterate.primal_values_[1] = 2.0;
	iterate.dual_eq_values_[0] = 2.78;
	iterate.dual_ieq_values_[0] = 4.66;
	
	iSQOIterate final(2,1,1);
	final.penalty_parameter_ = 1.0e-1;
	final.primal_values_[0] = .822875656;
	final.primal_values_[1] = .911437828;
	final.dual_eq_values_[0] = final.penalty_parameter_*1.846589027861980e+00;
	final.dual_ieq_values_[0] = final.penalty_parameter_*1.594493103554523e+00;
	
	     
    
	PenaltyFunction penfunc(&problem);
	
	cout << "current penalty is: " << penfunc(sillystart) << endl;
	cout << "current penalty is: " << penfunc(iterate) << endl;
	cout << "current penalty is: " << penfunc(final) << endl;
	
	NlpResidualFunction residual_func(&problem);
	cout << "residual @ final: " << residual_func(final) << endl;
	
	
	
	return 0;
}