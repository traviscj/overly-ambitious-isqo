
#include <vector>

#include "iterate.hh"
#include "nlp.hh"
#include "nlp_hs014.hh"

Hs014::Hs014() : Nlp(2,1,1) {
}

iSQOIterate Hs014::initial() {
	iSQOIterate tmp(2,1,1, 1e-1);
    std::vector<double> primal_values;
    primal_values.push_back(2.0);
    primal_values.push_back(1.0);
    std::vector<double> dual_eq_values;
	dual_eq_values.push_back(-1.0);
    std::vector<double> dual_ieq_values;
	dual_ieq_values.push_back(1.0);
    
    tmp.assign_primal(primal_values);
    tmp.assign_dual_eq(dual_eq_values);
    tmp.assign_dual_ieq(dual_ieq_values);
	return tmp;
}
double Hs014::objective(const iSQOIterate &iterate) {
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	return (primal_values[0] - 2)*(primal_values[0] - 2) + (primal_values[1]-1)*(primal_values[1]-1);
}
std::vector<double> Hs014::constraints_equality(const iSQOIterate &iterate) {
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
    
	std::vector<double> retval(num_dual_eq());
	retval[0] = primal_values[0] - 2*primal_values[1] + 1;
	return retval;
}
std::vector<double> Hs014::constraints_inequality(const iSQOIterate &iterate) {
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	std::vector<double> retval(num_dual_ieq());
	retval[0] = .25*primal_values[0]*primal_values[0] + primal_values[1]*primal_values[1] - 1;
	return retval;
}

std::vector<double> Hs014::objective_gradient(const iSQOIterate &iterate) {
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	std::vector<double> ret_gradient(num_primal());
	ret_gradient[0] = 2*(primal_values[0] - 2);
	ret_gradient[1] = 2*(primal_values[1] - 1);
	return ret_gradient;
}
std::shared_ptr<matrix_base_class> Hs014::constraints_equality_jacobian(const iSQOIterate &iterate){
	std::shared_ptr<dense_matrix> ret_jacobian(new dense_matrix(num_dual_eq(), num_primal()));
	ret_jacobian->set(0,0, 1.0);
	ret_jacobian->set(0,1,-2.0);
	return ret_jacobian;
}
std::shared_ptr<matrix_base_class> Hs014::constraints_inequality_jacobian(const iSQOIterate &iterate){
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	std::shared_ptr<dense_matrix> ret_jacobian(new dense_matrix(num_dual_ieq(), num_primal()));
	ret_jacobian->set(0,0,0.5*primal_values[0]);
	ret_jacobian->set(0,1,2.0*primal_values[1]);
	return ret_jacobian;
}

std::shared_ptr<matrix_base_class> Hs014::lagrangian_hessian(const iSQOIterate &iterate) {
	std::vector<double> primal_values(num_primal());
    primal_values.assign(iterate.get_primal_values()->begin(), iterate.get_primal_values()->end());
	std::vector<double> dual_ieq_values(num_primal());
    primal_values.assign(iterate.get_dual_ieq_values()->begin(), iterate.get_dual_ieq_values()->end());
    
	std::shared_ptr<dense_matrix> hess(new dense_matrix(num_primal(), num_primal()));
	if (iterate.get_penalty_parameter() > 0) {
		hess->set(0,0,iterate.get_penalty_parameter()*2.0 + 0.5*dual_ieq_values[0]); //   /iterate.penalty_parameter_
		hess->set(1,1,iterate.get_penalty_parameter()*2.0 + 2.0*dual_ieq_values[0]);	//	/iterate.penalty_parameter_
	} else {
		hess->set(0,0,0.5*dual_ieq_values[0]);
		hess->set(1,1,2.0*dual_ieq_values[1]);
	}
	return hess;
}
