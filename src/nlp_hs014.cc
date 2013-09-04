
#include <vector>

#include "iterate.hh"
#include "nlp.hh"
#include "nlp_hs014.hh"

Hs014::Hs014() : Nlp(2,1,1) {
}

iSQOIterate Hs014::initial() {
	iSQOIterate tmp(2,1,1);
	tmp.primal_values_[0] = 2.0;
	tmp.primal_values_[1] = 1.0;
	tmp.dual_eq_values_[0] = -1.0;
	tmp.dual_ieq_values_[0] = 1.0;
	return tmp;
}
double Hs014::objective(const iSQOIterate &iterate) {
	return (iterate.primal_values_[0] - 2)*(iterate.primal_values_[0] - 2) + (iterate.primal_values_[1]-1)*(iterate.primal_values_[1]-1);
}
std::vector<double> Hs014::constraints_equality(const iSQOIterate &iterate) {
	std::vector<double> retval(iterate.num_dual_eq_);
	retval[0] = iterate.primal_values_[0] - 2*iterate.primal_values_[1] + 1;
	return retval;
}
std::vector<double> Hs014::constraints_inequality(const iSQOIterate &iterate) {
	std::vector<double> retval(iterate.num_dual_ieq_);
	retval[0] = .25*iterate.primal_values_[0]*iterate.primal_values_[0] + iterate.primal_values_[1]*iterate.primal_values_[1] - 1;
	return retval;
}

std::vector<double> Hs014::objective_gradient(const iSQOIterate &iterate) {
	std::vector<double> ret_gradient(iterate.num_primal_);
	ret_gradient[0] = 2*(iterate.primal_values_[0] - 2);
	ret_gradient[1] = 2*(iterate.primal_values_[1] - 1);
	return ret_gradient;
}
std::shared_ptr<matrix_base_class> Hs014::constraints_equality_jacobian(const iSQOIterate &iterate){
	std::shared_ptr<dense_matrix> ret_jacobian(new dense_matrix(iterate.num_dual_ieq_, iterate.num_primal_));
	ret_jacobian->set(0,0, 1.0);
	ret_jacobian->set(0,1,-2.0);
	return ret_jacobian;
}
std::shared_ptr<matrix_base_class> Hs014::constraints_inequality_jacobian(const iSQOIterate &iterate){
	std::shared_ptr<dense_matrix> ret_jacobian(new dense_matrix(iterate.num_dual_ieq_, iterate.num_primal_));
	ret_jacobian->set(0,0,0.5*iterate.primal_values_[0]);
	ret_jacobian->set(0,1,2.0*iterate.primal_values_[1]);
	return ret_jacobian;
}

std::shared_ptr<matrix_base_class> Hs014::lagrangian_hessian(const iSQOIterate &iterate) {
	std::shared_ptr<dense_matrix> hess(new dense_matrix(iterate.num_primal_, iterate.num_primal_));
	if (iterate.penalty_parameter_ > 0) {
		hess->set(0,0,iterate.penalty_parameter_*2.0 + 0.5*iterate.dual_ieq_values_[0]); //   /iterate.penalty_parameter_
		hess->set(1,1,iterate.penalty_parameter_*2.0 + 2.0*iterate.dual_ieq_values_[0]);	//	/iterate.penalty_parameter_
	} else {
		hess->set(0,0,0.5*iterate.dual_ieq_values_[0]);
		hess->set(1,1,2.0*iterate.dual_ieq_values_[1]);
	}
	return hess;
}
