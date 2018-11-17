/**
    Binary header
    Binary.h
    Purpose: Functions calculate binary arrays and strings

    @author J.M.G. Salmerón
    @version 1.0 21/09/2018
    
    Created on: 21/09/2018
*/

#ifndef BINARY_H_
#define BINARY_H_

#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <bitset> // int_to_bin

#include "Types.h"
/*
class Binary {
    std::map<int, std::vector<std::string>> binaries;
    
public:
    Binary(int n);
    std::vector<std::string> k_ones_in_n_bits(int k, int n);
};
*/
std::vector<ULLI> gen_next_level(std::vector<ULLI> level, int n);
int count1s(std::string z);
std::vector<std::string> get_nextzs(std::string z);
void check_subfacets_live(std::string z, std::vector<ULLI> &next_level, std::vector<bool> &level_checked);
void check_nonsubfacets_live(std::string z, std::string new_z, std::vector<ULLI> &next_level, std::vector<bool> &next_level_checked);
ULLI bin_to_int(std::string z);
std::string int_to_bin(ULLI k, int n);

#endif /* BINARY_H_ */