/** @file
 *****************************************************************************
 How many additions in one second?
 *****************************************************************************
 * @author     Johannes Reinhart, ILS
 *****************************************************************************/

#include <iostream>
#include <math.h>

#include "depends/libff/libff/common/profiling.hpp"
#include "depends/libff/libff/common/utils.hpp"
#include "depends/libff/libff/common/default_types/ec_pp.hpp"

typedef libff::G1<libff::default_ec_pp> G1;
typedef libff::G2<libff::default_ec_pp> G2;
typedef libff::Fr<libff::default_ec_pp> Fr;
typedef libff::Fq<libff::default_ec_pp> Fq;
typedef libff::GT<libff::default_ec_pp> GT;

void profile_g1_additions(size_t iterations){
    G1 a = G1::random_element();
    G1 result1 = G1::G1_zero;
    G1 result2 = Fr(iterations)*a;
    libff::enter_block("G1 additions");

    for(size_t i = 0; i < iterations; i++){
        result1 = result1 + a;
    }
    libff::leave_block("G1 additions");
    bool success = (result1 == result2);
    std::cout << "Check: " << success << " result1: " << result1 << " result2: " << result2 << std::endl;
    a.print_coordinates();
    result1.print_coordinates();
    result2.print_coordinates();
}

void profile_g2_additions(size_t iterations){
    G2 a = G2::random_element();
    G2 result1 = G2::G2_zero;
    G2 result2 = Fr(iterations)*a;
    libff::enter_block("G2 additions");

    for(size_t i = 0; i < iterations; i++){
        result1 = result1 + a;
    }
    libff::leave_block("G2 additions");
    bool success = (result1 == result2);
    std::cout << "Check: " << success << " result1: " << result1 << " result2: " << result2 << std::endl;
    a.print_coordinates();
    result1.print_coordinates();
    result2.print_coordinates();
}

void profile_gt_additions(size_t iterations){
    GT a = GT::random_element();
    GT result1 = GT::zero();
    GT result2 = Fq(iterations)*a;
    libff::enter_block("GT additions");

    for(size_t i = 0; i < iterations; i++){
        result1 = result1 + a;
    }
    libff::leave_block("GT additions");
    bool success = (result1 == result2);
    std::cout << "Check: " << success << " result1: " << result1 << " result2: " << result2 << std::endl;
    a.print();
    result1.print();
    result2.print();
}


int main(int argc, char **argv) {
    libff::default_ec_pp::init_public_params();
    libff::start_profiling();

    size_t iterations;
    if (argc != 2){
        std::cout << "Usage: log10(num_iterations)" << std::endl;
        exit(-1);
    }

    iterations = std::pow(10, atoi(argv[1]));
    std::cout << "Profiling point additions: " << iterations << std::endl;

    profile_g1_additions(iterations);
    profile_g2_additions(iterations);
    profile_gt_additions(iterations);
}
