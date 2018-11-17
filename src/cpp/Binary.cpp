/**
    Binary body
    Binary.cpp
    Purpose: Functions calculate binary arrays and strings

    @author J.M.G. Salmerón
    @version 1.0 21/09/2018
    
    Created on: 21/09/2018
*/

#include "Binary.h"
/*
Binary::Binary(int n) {
	// TODO: Calculate k over n to know the size of a std::string* array
	std::string str;

    // construct N-digit binary number filled with all 0's
    int j = n;
    while (j--)
        str.push_back('0');

    // print all numbers with k-bit set together in ascending order
    for (int k = 1; k <= n; k++) {
        // set last k bits to '1'
        str[n - k] = '1';
        std::string curr = str;
        
        // use std::next_permutation to print string lexicographically
        do {
			binaries[k].push_back(curr);
        } while (next_permutation(curr.begin(), curr.end()));
    }
}
*/
/**
 * Get binary numbers of k 1s in n bits
 * 
 * @param k: Number of ones in the string
 * @param n: Number of bits of the string
 * @return Vector of strings
 */
/*
std::vector<std::string> Binary::k_ones_in_n_bits(int k, int n) {
	return binaries[k];
}
*/
/**
 * Return next level of binary strings from actual level
 * 
 * @param level: Actual level
 * @return Vector of strings
 */
std::vector<ULLI> gen_next_level(std::vector<ULLI> level, int n) {
	std::vector<ULLI> next_level;
    std::string z, new_z;

    for (unsigned int i = 0; i < size(level); i++) {
        z = int_to_bin(level[i], n);
        for (unsigned int j = 0; j < size(z); j++) {
            if (z[j] == '1') {
                new_z = z;
                new_z[j] = '0';
                next_level.push_back(bin_to_int(new_z));
            }
        }
    }

    sort(next_level.begin(), next_level.end());
    next_level.erase(unique(next_level.begin(), next_level.end()), next_level.end());

    return next_level;
}

/**
 * Count 1s in string
 * 
 * @param z: Bit string
 * @return Number of 1s
 */
int count1s(std::string z) {
    return std::count(z.begin(), z.end(), '1');
}

/**
 * Get subfacets of this facet
 * 
 * @return Vector of bit string with subfacets
 */
std::vector<std::string> get_nextzs(std::string z) {
    std::vector<std::string> sub_zs;

    for (unsigned int i = 0; i < size(z); i++) {
        if (z[i] == '1') {
            std::string sub_z(z);
            sub_z[i] = '0';
            sub_zs.push_back(sub_z);
        }
    }

    return sub_zs;
}

/**
 * Check subfacets in level who are children of F_k
 * 
 * @param z: Facet to check sons
 * @param next_level: Vector on bit strings
 * @param next_level_checked: Vector of checked facets
 */
void check_subfacets_live(std::string z, std::vector<ULLI> &next_level, std::vector<bool> &next_level_checked) {
    std::vector<std::string> zs = get_nextzs(z);

    for (unsigned int nl_i = 0; nl_i < size(next_level); nl_i++)
        for (unsigned int zs_i = 0; zs_i < size(zs); zs_i++)
            if (bin_to_int(zs[zs_i]) == next_level[nl_i])
                next_level_checked[nl_i] = true;
}

/**
 * Check subfacets in level who are children of F_k and are not new_z
 * 
 * Significa que he reducido la dimension, el resto de facetas de Fk las quito
 * Marcar las sub-facetas que no tengan x_i = 0 como checkeadas
 * 
 * @param z: Facet to check sons
 * @param new_z: Facet to check sons
 * @param next_level: Vector on bit strings
 * @param next_level_checked: Vector of checked facets
 */
void check_nonsubfacets_live(std::string z, std::string new_z, std::vector<ULLI> &next_level, std::vector<bool> &next_level_checked) {
    std::vector<std::string> zs = get_nextzs(z);

    for (unsigned int nl_i = 0; nl_i < size(next_level); nl_i++)
        for (unsigned int zs_i = 0; zs_i < size(zs); zs_i++)
            if (bin_to_int(zs[zs_i]) == next_level[nl_i] && bin_to_int(new_z) != next_level[nl_i])
                next_level_checked[nl_i] = true;
}

/**
 * Binary number (string) to integer
 * 
 * @param z: Binary string
 * @return Int number
 */
ULLI bin_to_int(std::string z) {
    return std::stoll(z, nullptr, 2); // C++11
}

/**
 * Integer number to binary string
 * 
 * @param k: Integer number
 * @return Binary string
 */
std::string int_to_bin(ULLI k, int n) {
    return std::bitset<64>(k).to_string().substr(64 - n, 64);
}
