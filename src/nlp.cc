
#include "nlp.hh"

Nlp::Nlp(int num_primal, int num_dual_eq, int num_dual_ieq) : num_primal_(num_primal), num_dual_eq_(num_dual_eq), num_dual_ieq_(num_dual_ieq) {
	// cout << "-- initializing Nlp" << endl; 
}
int Nlp::num_primal() { return num_primal_; }
int Nlp::num_dual() {return num_dual_eq_+num_dual_ieq_; }
int Nlp::num_dual_eq() { return num_dual_eq_; }
int Nlp::num_dual_ieq() { return num_dual_ieq_; }
