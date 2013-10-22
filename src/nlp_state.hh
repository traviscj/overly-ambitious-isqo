// Copyright (C) 2013 Travis C. Johnson (traviscj@traviscj.com)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author:  Travis C. Johnson (traviscj@traviscj.com)

#ifndef __GUARD_NLP_STATE_HH
#define __GUARD_NLP_STATE_HH

#include "nlp.hh"

//! \brief a class for holding some Nlp state. derived classes override operator()
// \todo implement some caching here.
class FunctionWithNLPState {
public:
	FunctionWithNLPState(Nlp &nlp) {
		nlp_ = &nlp;
	}
protected:
	Nlp *nlp_;
};

#endif
