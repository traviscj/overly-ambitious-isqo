
#include <iostream>
#include <map>

#include "isqo_config.hh"

// TODO should be a singleton??
iSQOControlPanel::iSQOControlPanel() {
    hessian_shift_print_level = HSPL_NONE;
    initial_penalty_parameter = 1e-1;
    qp_print_level=QPL_NONE;
    termination_tolerance=1e-6;
    
    // setup the maps from config strings -> *_t
    map_hessian_shift_print_levels["none"] = HSPL_NONE;
    map_hessian_shift_print_levels["full"] = HSPL_FULL;
    map_qp_print_levels["none"] = QPL_NONE;
    map_qp_print_levels["low"] = QPL_LOW;
    map_qp_print_levels["medium"] = QPL_MEDIUM;
    map_qp_print_levels["high"] = QPL_HIGH;
    map_qp_print_levels["table"] = QPL_TABLE;
} 

template <typename T, typename S>
void iSQOControlPanel::update_junk(std::istream &is, std::map<std::string, T> update_map, S &value_to_update) {
    std::string scratch;
    is >> scratch;
    if (update_map.find(scratch) != update_map.end()) {
        value_to_update = update_map[scratch];
    } else {
        std::cerr << "unknown " << "??? option specified: [" << scratch << "]."<< std::endl;
    }
}
void iSQOControlPanel::update_settings(std::istream &is) {
    std::string identifier;
    // string scratch;
    while (!is.eof()) {
        is >> identifier;
        if (identifier.compare("hessian_shift_print_level:") == 0) {
            update_junk(is, map_hessian_shift_print_levels, hessian_shift_print_level);
        } else if (identifier.compare("initial_penalty_parameter:") == 0) {
            is >> initial_penalty_parameter;
        } else if (identifier.compare("qp_print_level:") == 0) {
            update_junk(is, map_qp_print_levels, qp_print_level);
        } else if (identifier.compare("termination_tolerance:") == 0) {
            is >> termination_tolerance;
        } else {
            std::cerr << "undefined identifier in input file: [" << identifier << "]." << std::endl;
        }
    }
}

void iSQOControlPanel::print() {
    std::cout << "All settings:" << std::endl
        << "hessian_shift_print_level_t hessian_shift_print_level = " << hessian_shift_print_level << std::endl
        << "qp_print_level_t qp_print_level = " << qp_print_level << std::endl
        << "double initial_penalty_parameter = " << initial_penalty_parameter << std::endl
        << "double termination_tolerance = " << termination_tolerance << std::endl
        << std::endl;
}
hessian_shift_print_level_t iSQOControlPanel::get_hessian_shift_print_level() { return hessian_shift_print_level; }
double iSQOControlPanel::get_initial_penalty_parameter() { return initial_penalty_parameter; }
qp_print_level_t iSQOControlPanel::get_qp_print_level() { return qp_print_level; }
double iSQOControlPanel::get_termination_tolerance() { return termination_tolerance; }

