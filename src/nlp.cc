
#include "nlp.hh"

Nlp::Nlp(size_t num_primal, size_t num_dual_eq, size_t num_dual_ieq) : num_primal_(num_primal), num_dual_eq_(num_dual_eq), num_dual_ieq_(num_dual_ieq) {
	// cout << "-- initializing Nlp" << endl; 
}

size_t Nlp::num_primal() { return num_primal_; }
size_t Nlp::num_dual() {return num_dual_eq_+num_dual_ieq_; }
size_t Nlp::num_dual_eq() { return num_dual_eq_; }
size_t Nlp::num_dual_ieq() { return num_dual_ieq_; }
