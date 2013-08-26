
#ifndef __GUARD_NLP_STATE_HH
#define __GUARD_NLP_STATE_HH

#include "nlp.hh"

class FunctionWithNLPState {
public:
	FunctionWithNLPState(Nlp &nlp) {
		nlp_ = &nlp;
	}
protected:
	Nlp *nlp_;
};

#endif
