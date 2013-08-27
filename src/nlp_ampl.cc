
#include <iostream>
#include <cassert>

#include "nlp_ampl.hh"
#include "utilities.hh"


// AmplNlp(char *stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
AmplNlp::AmplNlp(std::string stub_str) : Nlp(-1,-1,-1), PRINT_(false) {
	if (PRINT_) std::cout << "Constructing an AmplNlp" << std::endl;
	ASL *asl;
	asl = ASL_alloc(ASL_read_pfgh);
	nerror_ = new int;
	*nerror_ = 0;
	// char * punt = stub_str.data();
	// if (PRINT_) 
		std::cout << "Filename: " << stub_str << std::endl;
	nl_ = jac0dim(&stub_str[0], (stub_str).length());
	if (PRINT_) 
		std::cout << "ran jac0dim: " << nl_ << std::endl;
		
	num_primal_ = n_var;
	X0_.resize(num_primal());
	X0 = &X0_[0];
	// A_vals = new double[nzc];
	// A_rownos = new int[nzc];
	// A_colstarts = new int[nzc];
	// ideally, we would use these....

	int status = pfgh_read(nl_, ASL_return_read_err | ASL_findgroups);

	if (PRINT_) std::cout << "X0: " << X0[0] << ", " << X0[1] << std::endl;
	
	for (size_t ampl_constraint_index=0; ampl_constraint_index<n_con; ++ampl_constraint_index) {
		if (PRINT_) std::cout << "On constraint ampl_constraint_index=" << ampl_constraint_index << ": ";
		bool equality_con = LUrhs[2*ampl_constraint_index]==LUrhs[2*ampl_constraint_index+1] && LUrhs[2*ampl_constraint_index]!=INFINITY;
	
		if (equality_con) {
			if (PRINT_) std::cout << "equality";
			if (PRINT_) std::cout << "(push eq)";
			equality_constraints_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] == -INFINITY && LUrhs[2*ampl_constraint_index+1] > -INFINITY) {
			if (PRINT_) std::cout << "upper bound";
			if (PRINT_) std::cout << "(push ineq)";
			inequality_constraints_upper_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY) {
			if (PRINT_) std::cout << "double-sided";
			if (PRINT_) std::cout << "(push ineq)";
			if (PRINT_) std::cout << "(push -ineq)";
			inequality_constraints_upper_.push_back(ampl_constraint_index);
			inequality_constraints_lower_.push_back(ampl_constraint_index);
		} else if (LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY) {
			if (PRINT_) std::cout << "lower bound only";
			if (PRINT_) std::cout << "(push -ineq)";
			inequality_constraints_lower_.push_back(ampl_constraint_index);
		} else {
			std::cout << "unsupported";
		}
		if (PRINT_) std::cout << std::endl;
	}
	
	
	
	std::vector<double> x(num_primal());
	// double *x = new double[num_primal()];
	for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
		x[primal_index] = X0_[primal_index];
	
	bool should_abort = false;
	for (size_t ampl_variable_index =0 ; ampl_variable_index < n_var; ++ampl_variable_index) {
		if (PRINT_) 
			std::cout << "ampl_variable_index=" << ampl_variable_index << ": l=" << LUv[2*ampl_variable_index] << " <= x_ampl_variable_index= " << x[ampl_variable_index] << " <= u=" << LUv[2*ampl_variable_index+1] << std::endl;
		bool lower_finite = (LUv[2*ampl_variable_index] != -INFINITY);
		bool upper_finite = (LUv[2*ampl_variable_index+1] != INFINITY);
		bool both_finite = lower_finite && upper_finite;
		bool equal = (LUv[2*ampl_variable_index] == LUv[2*ampl_variable_index+1]);
		if (both_finite && equal) {
			variable_equality_.push_back(ampl_variable_index);
		} else if (both_finite) {
			variable_bound_lower_.push_back(ampl_variable_index);
			variable_bound_upper_.push_back(ampl_variable_index);
		} else if (lower_finite) {
			variable_bound_lower_.push_back(ampl_variable_index);
		} else if (upper_finite) {
			variable_bound_upper_.push_back(ampl_variable_index);
		}			
	}
	
    PRINT_ = true;
	std::vector<double> con(n_con);
	conval(&x[0], &con[0], nerror_);
	for (size_t ampl_constraint_index =0 ; ampl_constraint_index < n_con; ++ampl_constraint_index) {
		if (PRINT_) std::cout << "ampl_constraint_index=" << ampl_constraint_index << ": l=" << LUrhs[2*ampl_constraint_index] << " <= c(x_k)= " << con[ampl_constraint_index] << " <= u=" << LUrhs[2*ampl_constraint_index+1] << std::endl;
	}
    PRINT_ = false;
	
	num_dual_eq_ = equality_constraints_.size();
	num_dual_ieq_ = inequality_constraints_lower_.size() + inequality_constraints_upper_.size() + variable_bound_lower_.size() + variable_bound_upper_.size();
	asl_ = asl;
}

AmplNlp::~AmplNlp() {
	ASL *asl = asl_;
	
	delete nerror_;
	ASL_free(&asl);
	fclose(nl_);
}

iSQOIterate AmplNlp::initial() {
	iSQOIterate initial(num_primal(), num_dual_eq(), num_dual_ieq());
	for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
		initial.primal_values_[primal_index] = X0_[primal_index];
		if (PRINT_) std::cout << "initial_" << primal_index << ": " << X0_[primal_index] <<std::endl;
	}
	std::vector<double> eq_values = constraints_equality(initial);
	for (size_t dual_eq_index=0; dual_eq_index < num_dual_eq(); ++dual_eq_index) {
		if (eq_values[dual_eq_index] < 0) {
			initial.dual_eq_values_[dual_eq_index] = -1.0;
		} else if (eq_values[dual_eq_index] > 0) {
			initial.dual_eq_values_[dual_eq_index] = +1.0;
		} else {
			initial.dual_eq_values_[dual_eq_index] = 0.0;
		}
	}

	std::vector<double> ieq_values = constraints_inequality(initial);
	for (size_t dual_ieq_index=0; dual_ieq_index < num_dual_ieq(); ++dual_ieq_index) {
		if (ieq_values[dual_ieq_index] < 0) {
			initial.dual_ieq_values_[dual_ieq_index] = 0.0;
		} else if (ieq_values[dual_ieq_index] > 0) {
			initial.dual_ieq_values_[dual_ieq_index] =+1.0;
		} else {
			initial.dual_ieq_values_[dual_ieq_index] = 0.0;
		}
	}
	return initial;
}
// zeroth order NLP quantities:
double AmplNlp::objective(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	
	real obj = objval(0, &x[0], nerror_);
	if (PRINT_) std::cout << "objective value(nerror = " << *nerror_ << "): " << obj << std::endl;
	return obj;
}

std::vector<double> AmplNlp::constraints_equality(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	// std::vector<double> con(num_dual());
	std::vector<double> con(n_con);
	conval(&x[0], &con[0], nerror_);
	if (PRINT_) std::cout << "equality_constraints:" << std::endl;
	std::vector<double> equality_constraint_evaluation(equality_constraints_.size());
	for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
		equality_constraint_evaluation[isqo_eq_constraint_index] = con[equality_constraints_[isqo_eq_constraint_index]] - LUrhs[2*equality_constraints_[isqo_eq_constraint_index]];
	}
	if (PRINT_) std::cout << "eq eval: " << equality_constraint_evaluation;
	return equality_constraint_evaluation;
}
std::vector<double> AmplNlp::constraints_inequality(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	
	std::vector<double> con(n_con);
	
	conval(&x[0], &con[0], nerror_);
	if (PRINT_) std::cout << "inequality_constraints:" << std::endl;
	
	std::vector<double> inequality_constraint_evaluation(num_dual_ieq());
	
	int isqo_ineq_constraint_index = 0;
	for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
		if (PRINT_) std::cout << "isqo ineq " << isqo_ineq_constraint_index << "; " << inequality_constraints_lower_[isqo_ineq_lower_constraint_index] << std::endl;
		size_t ampl_constraint_index = inequality_constraints_lower_[isqo_ineq_lower_constraint_index];
		size_t ampl_constraint_value_lower = 2*ampl_constraint_index;
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = -con[ampl_constraint_index] + LUrhs[ampl_constraint_value_lower];
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
		
		size_t ampl_constraint_index = inequality_constraints_upper_[isqo_ineq_upper_constraint_index];
		size_t	ampl_constraint_value_upper = 2*ampl_constraint_index+1;
		if (PRINT_) std::cout << "isqo ineq " << isqo_ineq_upper_constraint_index << " is AMPL constraint " << ampl_constraint_index << std::endl;
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = con[ampl_constraint_index] - LUrhs[ampl_constraint_value_upper];
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]];
		// std::cout 
		// 	<< "variable lower bound " << isqo_ineq_lower_variable_index 
		// 	<< "; ampl bound: " <<  variable_bound_lower_[isqo_ineq_lower_variable_index]
		// 	<< "; with value: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]]
		// 	<< "; @x_k: " << LUv[2*variable_bound_lower_[isqo_ineq_lower_variable_index]] - x[variable_bound_lower_[isqo_ineq_lower_variable_index]]
		// 	<< std::endl;
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
		// std::cout << "variable upper bound " << isqo_ineq_upper_variable_index << "" << std::endl;
		inequality_constraint_evaluation[isqo_ineq_constraint_index] = -LUv[2*variable_bound_lower_[isqo_ineq_upper_variable_index]+1] + x[variable_bound_lower_[isqo_ineq_upper_variable_index]];
		++isqo_ineq_constraint_index;
	}
	assert(isqo_ineq_constraint_index == num_dual_ieq());
	
	
	if (PRINT_) std::cout <<"ineq eval" << inequality_constraint_evaluation;
	return inequality_constraint_evaluation;
}

// first order NLP quantities:
std::vector<double> AmplNlp::objective_gradient(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> return_gradient(num_primal());
	
	std::vector<double> x(num_primal());
	for (size_t primal_index=0; primal_index < num_primal(); ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	
	objgrd(0, &x[0], &return_gradient[0], nerror_);
	if (PRINT_) std::cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[0] << std::endl;
	if (PRINT_) std::cout << "objective gradient(nerror = " << *nerror_ << "): " << return_gradient[1] << std::endl;
	
	return return_gradient;
}
matrix AmplNlp::constraints_equality_jacobian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	
	matrix equality_constraint_jacobian(equality_constraints_.size(), n_var);
	if (PRINT_) std::cout << "equal: " << equality_constraints_ << std::endl;
	std::vector<double> G(num_primal());
	for (size_t isqo_eq_constraint_index=0; isqo_eq_constraint_index < equality_constraints_.size(); ++isqo_eq_constraint_index) {
		congrd(equality_constraints_[isqo_eq_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			equality_constraint_jacobian.set(isqo_eq_constraint_index, primal_index, G[primal_index]);
		}
		
	}
	if (PRINT_) std::cout << "equality jacobian: [" << std::endl;
	if (PRINT_) std::cout << equality_constraint_jacobian;
	if (PRINT_) std::cout << "]" << std::endl;

	return equality_constraint_jacobian;
}
qpOASES::SparseMatrix AmplNlp::constraints_equality_jacobian_sparse(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index)
		x[primal_index] = iterate.primal_values_[primal_index];
	
    asl->i.congrd_mode = 0;
    double J[nzc];
    jacval(&x[0], &J[0], nerror_);
    std::cout << "J: " << std::endl;
    for (size_t current_nz=0; current_nz < nzc; ++current_nz) {
        std::cout << "current_nz=" << current_nz << ", val=" << J[current_nz] << std::endl;
    }
    
    cgrad *cg;
    size_t my_count_nz =0;
    
    std::cout << "equ: " << equality_constraints_ << std::endl;
    std::cout << "ineq lower: " << inequality_constraints_lower_ << std::endl;
    std::cout << "ineq upper: " << inequality_constraints_upper_ << std::endl;
    
    std::vector<int> row_coordinate(nzc);
    std::vector<int> col_coordinate(nzc);
    std::vector<double> value(nzc);

    std::cout << "FULL TABLE: " << std::endl;
    // the number of nonzero entries in each type:
    size_t num_eq_nnz=0, num_ieq_lower_nnz=0, num_ieq_upper_nnz=0;
    for (size_t ampl_constraint_index=0; ampl_constraint_index<n_con; ++ampl_constraint_index) {
        for (cg = Cgrad[ampl_constraint_index]; cg; cg = cg->next){
            char constraint_type;
            if ( LUrhs[2*ampl_constraint_index] == LUrhs[2*ampl_constraint_index+1] ) {
                constraint_type = 'E';
                ++num_eq_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY ) { 
                constraint_type = 'B';
                ++num_ieq_lower_nnz;
                ++num_ieq_upper_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] == -INFINITY && LUrhs[2*ampl_constraint_index+1] < INFINITY ) { 
                constraint_type = 'U';
                ++num_ieq_upper_nnz;
            } else if ( LUrhs[2*ampl_constraint_index] > -INFINITY && LUrhs[2*ampl_constraint_index+1] == INFINITY ) { 
                constraint_type = 'L';
                ++num_ieq_lower_nnz;
            } else {
                assert(false);
            }
            
            printf("index: %3d: conindex: %3d[%9.2e <= c_i <= %9.2e ==> %c]: varno: %3d, coef: %9.2e, goff: %3d, val: %9.2e\n", 
                    my_count_nz, ampl_constraint_index, LUrhs[2*ampl_constraint_index], LUrhs[2*ampl_constraint_index+1], constraint_type,
                    cg->varno, cg->coef, cg->goff, J[cg->goff]);
            row_coordinate[cg->goff] = ampl_constraint_index;
            col_coordinate[cg->goff] = cg->varno;
            value[cg->goff] = J[cg->goff];
            // ++col_lengths_eq[cg->varno];
            // ++current_eq_index;
            ++my_count_nz;
        }
    }
    for (size_t ampl_variable_index=0; ampl_variable_index < n_var; ++ampl_variable_index) {
        char variable_type;
        if ( LUv[2*ampl_variable_index] == LUv[2*ampl_variable_index+1] ) {
            variable_type = 'E';
        } else if ( LUv[2*ampl_variable_index] > -INFINITY && LUv[2*ampl_variable_index+1] < INFINITY ) {
            variable_type = 'B';
        } else if ( LUv[2*ampl_variable_index] == -INFINITY && LUv[2*ampl_variable_index+1] < INFINITY ) {
            variable_type = 'U';
        } else if ( LUv[2*ampl_variable_index] > -INFINITY && LUv[2*ampl_variable_index+1] == INFINITY ) {
            variable_type = 'L';
        } else {
            variable_type='F';
        }
        printf("variable index: %3d [%9.2e <= x_i <= %9.2e ==> %c]\n", ampl_variable_index, LUv[2*ampl_variable_index],LUv[2*ampl_variable_index+1],variable_type);
    }
    typedef enum {EQUALITY=1, LOWER=2, UPPER=4} constraint_type;
    // constraint_types holds one of the above enums depending on what type.
    std::vector<size_t> constraint_types(n_con);
    // constraint_index holds the index of the constraint
    std::vector<size_t> constraint_indices(n_con);
        
    for (size_t current_eq_index=0; current_eq_index < equality_constraints_.size(); ++current_eq_index) {
        constraint_types[equality_constraints_[current_eq_index]] |= EQUALITY;
        constraint_indices[equality_constraints_[current_eq_index]] = current_eq_index;
    }
    for (size_t current_ieq_lower_index=0; current_ieq_lower_index < inequality_constraints_lower_.size(); ++current_ieq_lower_index) {
        constraint_types[inequality_constraints_lower_[current_ieq_lower_index]] |= LOWER;
        constraint_indices[inequality_constraints_lower_[current_ieq_lower_index]] = current_ieq_lower_index;
    }
    for (size_t current_ieq_upper_index=0; current_ieq_upper_index < inequality_constraints_upper_.size(); ++current_ieq_upper_index) {
        constraint_types[inequality_constraints_upper_[current_ieq_upper_index]] |= UPPER;
        constraint_indices[inequality_constraints_upper_[current_ieq_upper_index]] = current_ieq_upper_index;
    }
        
    std::cout << "contype: " << constraint_types << std::endl;
    std::cout << "con ind: " << constraint_indices << std::endl;

    std::cout << "row_coordinate: " << row_coordinate << std::endl;
    std::cout << "col_coordinate: " << col_coordinate << std::endl;
    std::cout << "value: " << value << std::endl;
    
    std::cout << "ASDF: n_var: " << n_var << std::endl;
    std::cout << "ASDFGHJ: equality_constraints_.size() = " << equality_constraints_.size() << std::endl;
    std::cout << "ASDFGHJ: inequality_constraints_lower_.size() = " << inequality_constraints_lower_.size() << std::endl;
    std::cout << "ASDFGHJ: inequality_constraints_upper_.size() = " << inequality_constraints_upper_.size() << std::endl;
    
    sparse_matrix eq_jacobian(equality_constraints_.size(), n_var, num_eq_nnz); 
    sparse_matrix ieq_lower_jacobian(inequality_constraints_lower_.size(), n_var, num_ieq_lower_nnz);
    sparse_matrix ieq_upper_jacobian(inequality_constraints_upper_.size(), n_var, num_ieq_upper_nnz);
    
    std::vector<int> col_lengths_eq(n_var);
    std::vector<int> col_lengths_ieq_lower(n_var);
    std::vector<int> col_lengths_ieq_upper(n_var);

    int current_eq_index = 0, current_ieq_lower_index = 0, current_ieq_upper_index = 0;
    
    for (size_t current_nonzero=0; current_nonzero < nzc; ++current_nonzero) {
        std::cout << "current_nonzero: " << current_nonzero << "; " << row_coordinate[current_nonzero] << ", " << col_coordinate[current_nonzero] << "; value: " << value[current_nonzero];
        if (constraint_types[row_coordinate[current_nonzero]] == EQUALITY) {
            std::cout << " - eq";
            eq_jacobian.vals_[current_eq_index] = value[current_nonzero];
            eq_jacobian.row_indices_[current_eq_index] = constraint_indices[row_coordinate[current_nonzero]];
            ++col_lengths_eq[col_coordinate[current_nonzero]];
            ++current_eq_index;
        }
        if (constraint_types[row_coordinate[current_nonzero]] == LOWER) {
            std::cout << " - ieq lower";
            ieq_lower_jacobian.vals_[current_ieq_lower_index] = -value[current_nonzero];
            ieq_lower_jacobian.row_indices_[current_ieq_lower_index] = constraint_indices[row_coordinate[current_nonzero]];
            ++col_lengths_ieq_lower[col_coordinate[current_nonzero]];
            ++current_ieq_lower_index;
        }
        if (constraint_types[row_coordinate[current_nonzero]] == UPPER) {
            std::cout << " - ieq upper";
            ieq_upper_jacobian.vals_[current_ieq_upper_index] = value[current_nonzero];
            ieq_upper_jacobian.row_indices_[current_ieq_upper_index] = constraint_indices[row_coordinate[current_nonzero]];
            ++col_lengths_ieq_upper[col_coordinate[current_nonzero]];
            ++current_ieq_upper_index;
        }
        std::cout << ".done" << std::endl;
    }
    std::cout << "check: " << current_eq_index << ", " << current_ieq_lower_index << ", " << current_ieq_upper_index << std::endl;
 
    std::cout << "col_lengths: " << col_lengths_eq << std::endl;
    size_t current_entry_eq = 0, current_entry_ieq_lower = 0, current_entry_ieq_upper = 0;
    for (size_t variable_index=0; variable_index<col_lengths_eq.size(); ++variable_index) {
        // EQUALITIES
        current_entry_eq += col_lengths_eq[variable_index];
        eq_jacobian.col_starts_[variable_index+1] = current_entry_eq;
        
        // INEQUALITIES: LOWER
        current_entry_ieq_lower += col_lengths_ieq_lower[variable_index];
        ieq_lower_jacobian.col_starts_[variable_index+1] = current_entry_ieq_lower;
        
        // INEQUALITIES: UPPER
        current_entry_ieq_upper += col_lengths_ieq_upper[variable_index];
        ieq_upper_jacobian.col_starts_[variable_index+1] = current_entry_ieq_upper;
        
    }
    
    std::cout << "nzc: " << nzc << "; my: " << my_count_nz << std::endl;
    std::cout << "equality ------------" << std::endl;
    std::cout << eq_jacobian << std::endl;
    
    matrix equality_constraint_jacobian = constraints_equality_jacobian(iterate);
    matrix inequality_constraint_jacobian = constraints_inequality_jacobian(iterate);
    PRINT_ = true;
    
    if (PRINT_) std::cout << "equality jacobian: [" << std::endl;
    if (PRINT_) std::cout << equality_constraint_jacobian;
    if (PRINT_) std::cout << "]" << std::endl;
    

    std::cout << "constraint inequality (lower) ------------" << std::endl;
    std::cout << ieq_lower_jacobian << std::endl;

    std::cout << "constraint inequality (upper) ------------" << std::endl;
    std::cout << ieq_upper_jacobian << std::endl;
    
    std::cout << "constraint inequality ------------" << std::endl;
    sparse_matrix constraints_ieq = vertical(ieq_lower_jacobian, ieq_upper_jacobian);
    std::cout << constraints_ieq << std::endl;
    
    
    std::cout << "\n\n\n=========" << std::endl;
    
    
    std::cout << "=========\n\n\n" << std::endl;
    
    sparse_matrix lower_var(n_var, variable_bound_lower_, -1.0);
    sparse_matrix upper_var(n_var, variable_bound_upper_, +1.0);
    std::cout << "lower_var: " << lower_var << std::endl;
    std::cout << "upper_var: " << upper_var << std::endl;
    
    
    sparse_matrix variable_bounds_jacobian = vertical(lower_var, upper_var);
    std::cout << "variable_bounds_jacobian: " << variable_bounds_jacobian << std::endl;
    
    std::cout << "=========\n\n\n" << std::endl;
    
    
    sparse_matrix full_ieq = vertical(constraints_ieq, variable_bounds_jacobian);
    std::cout << "full inequality ------------" << std::endl;
    std::cout << full_ieq << std::endl;
    
    qpOASES::SparseMatrix testmat_ieq(full_ieq.num_rows(), n_var, &full_ieq.row_indices_[0], &full_ieq.col_starts_[0], &full_ieq.vals_[0]);
    double *testmat_ieq_full = testmat_ieq.full();
    for (int i=0; i<200; i++) std::cout << ", " << testmat_ieq_full[i];
    std::cout << std::endl;
    std::cout << "testmat_ieq.full(): " << std::endl;
    // qpOASES::DenseMatrix dense_testmat_ieq(full_ieq.num_rows(), n_var, )
    
    if (PRINT_) std::cout << "inequality jacobian: [" << std::endl;
    if (PRINT_) std::cout << inequality_constraint_jacobian;
    if (PRINT_) std::cout << "]" << std::endl;
    
    for (int i=0; i<200; i++) assert(testmat_ieq_full[i] == inequality_constraint_jacobian.data_[i]);
    
    qpOASES::SparseMatrix testmat(n_con, n_var, &eq_jacobian.row_indices_[0], &eq_jacobian.col_starts_[0], &eq_jacobian.vals_[0]);
    // qpOASES::DenseMatrix testmat_dense(testmat);
    // testmat_dense.print();
    
    PRINT_ = false;

	return testmat;
}
matrix AmplNlp::constraints_inequality_jacobian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> x(num_primal());
	for (size_t primal_index=0; primal_index < iterate.num_primal_; ++primal_index){
		x[primal_index] = iterate.primal_values_[primal_index];
	}
	
	matrix inequality_constraint_jacobian(num_dual_ieq(), n_var);
	size_t isqo_ineq_constraint_index=0;
	if (PRINT_) std::cout << "lower: " << inequality_constraints_lower_ << std::endl;
	if (PRINT_) std::cout << "upper: " << inequality_constraints_upper_ << std::endl;
	std::vector<double> G(num_primal());
	for (size_t isqo_ineq_lower_constraint_index=0; isqo_ineq_lower_constraint_index < inequality_constraints_lower_.size(); ++isqo_ineq_lower_constraint_index) {
		congrd(inequality_constraints_lower_[isqo_ineq_lower_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, -G[primal_index]);
		}
		++isqo_ineq_constraint_index;
	}
	for (size_t isqo_ineq_upper_constraint_index=0; isqo_ineq_upper_constraint_index < inequality_constraints_upper_.size(); ++isqo_ineq_upper_constraint_index) {
		congrd(inequality_constraints_upper_[isqo_ineq_upper_constraint_index], &x[0], &G[0], nerror_);
		for (size_t primal_index=0; primal_index < num_primal(); ++primal_index) {
			inequality_constraint_jacobian.set(isqo_ineq_constraint_index, primal_index, G[primal_index]);
		}
		++isqo_ineq_constraint_index;
	}
    // set up Jacobian for variable bounds.
    // LOWER VARIABLE BOUNDS:
	for (size_t isqo_ineq_lower_variable_index=0; isqo_ineq_lower_variable_index < variable_bound_lower_.size(); ++isqo_ineq_lower_variable_index) {
		inequality_constraint_jacobian.set(isqo_ineq_constraint_index, variable_bound_lower_[isqo_ineq_lower_variable_index], -1.0);
		++isqo_ineq_constraint_index;
	}
    // UPPER VARIABLE BOUNDS:
	for (size_t isqo_ineq_upper_variable_index=0; isqo_ineq_upper_variable_index < variable_bound_upper_.size(); ++isqo_ineq_upper_variable_index) {
		inequality_constraint_jacobian.set(isqo_ineq_constraint_index, variable_bound_upper_[isqo_ineq_upper_variable_index], +1.0);
		++isqo_ineq_constraint_index;
	}
    
	if (PRINT_) 
	{
		std::cout << "======" << std::endl;
		std::cout << inequality_constraint_jacobian << std::endl;
		std::cout << "======" << std::endl;
	}
			
	return inequality_constraint_jacobian;
}

// second order NLP quantities:
matrix AmplNlp::lagrangian_hessian(const iSQOIterate &iterate) {
	ASL *asl = asl_;
	std::vector<double> H(num_primal() * num_primal());
	std::vector<double> OW(1);
	OW[0] = iterate.penalty_parameter_;
	std::vector<double> Y(n_con);
	
	// funnel dual_{,i}eq_values_ into Y to pass to AMPL.
	for (size_t dual_eq_index=0; dual_eq_index < iterate.dual_eq_values_.size(); ++dual_eq_index) {
		Y[equality_constraints_[dual_eq_index]] = iterate.dual_eq_values_[dual_eq_index];// / iterate.penalty_parameter_;
	}
	
	size_t dual_ineq_value_index = 0;
	for (size_t dual_ineq_lower_index=0; dual_ineq_lower_index < inequality_constraints_lower_.size(); ++dual_ineq_lower_index) {
		Y[inequality_constraints_lower_[dual_ineq_lower_index]] = -iterate.dual_ieq_values_[dual_ineq_value_index];
		dual_ineq_value_index++;
	}
	for (size_t dual_ineq_upper_index=0; dual_ineq_upper_index < inequality_constraints_upper_.size(); ++dual_ineq_upper_index) {
		Y[inequality_constraints_upper_[dual_ineq_upper_index]] = iterate.dual_ieq_values_[dual_ineq_value_index];
		dual_ineq_value_index++;
	}
	if (PRINT_) std::cout << "Y[0] = " << Y[0] << std::endl;
	if (PRINT_) std::cout << "Y[1] = " << Y[1] << std::endl;
	if (PRINT_) std::cout << "Y[2] = " << Y[2] << std::endl;
	if (PRINT_) std::cout << "Y[3] = " << Y[3] << std::endl;		
	
	// In this call:
	//	 - H : OUTPUT : vector holding entries of full hessian.
	//	 - n_var : INPUT : # of variables or something about the stride, maybe.
	//	 - 0 : INPUT : index of the desired objective.
	//	 - OW : INPUT : multipliers for objective function ("objective weights")
	//	 - Y : INPUT : lagrange multipliers for constraints.
	fullhes(&H[0], n_var, 0, &OW[0], &Y[0]);
	if (iterate.penalty_parameter_ == 0) {
		std::vector<double> H_f_only(num_primal()*num_primal());
		std::vector<double> OW_f_only(1);
		OW_f_only[0] = 1.0;
		std::vector<double> Y_f_only(n_con);
		for (size_t dual_index=0; dual_index < n_con; ++dual_index) Y[dual_index] = 0.0;
		fullhes(&H_f_only[0], n_var, 0, &OW_f_only[0], &Y_f_only[0]);
		
		for (size_t row_index=0; row_index < n_var; ++row_index)
			for (size_t col_index=0; col_index < n_var; ++col_index) 
				H[row_index*n_var + col_index] -= H_f_only[row_index*n_var + col_index];
	}
	
	matrix return_hessian(n_var, n_var);
	for (size_t row_index=0; row_index < n_var; ++row_index) {
		for (size_t col_index=0; col_index < n_var; ++col_index) {
			return_hessian.set(row_index, col_index, H[row_index*n_var + col_index]);
		}
	}
	return return_hessian;
}


std::size_t AmplNlp::num_lower_ieq() {
	return inequality_constraints_lower_.size();
}
std::size_t AmplNlp::num_upper_ieq() {
	return inequality_constraints_upper_.size();
}

