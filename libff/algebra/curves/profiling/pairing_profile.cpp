/**
 *****************************************************************************
 Profile pairings.

 *****************************************************************************/
#include <cstdio>
#include <vector>

#include "depends/libff/libff/algebra/curves/bn254/bn254_pp.hpp"
#include "depends/libff/libff/algebra/curves/bn183/bn183_pp.hpp"
#include "depends/libff/libff/algebra/curves/bn124/bn124_pp.hpp"
#include "depends/libff/libff/algebra/curves/edwards181/edwards181_pp.hpp"
#include "depends/libff/libff/algebra/curves/edwards97/edwards97_pp.hpp"
#include "depends/libff/libff/algebra/curves/edwards58/edwards58_pp.hpp"
#include "depends/libff/libff/algebra/scalar_multiplication/multiexp.hpp"
#include "depends/libff/libff/common/rng.hpp"

using namespace libff;

using std::size_t;


template <typename T>
using test_instances_t = std::vector<std::pair<G1<T>, G2<T>>>;

template<typename pp>
test_instances_t<pp> generate_group_elements(size_t count)
{
    test_instances_t<pp> result(count);

    for (size_t i = 0; i < count; i++) {
        G1<pp> g1 = G1<pp>::random_element();
        G2<pp> g2 = G2<pp>::random_element();

        result[i].first = g1;
        result[i].second = g2;
    }

    return result;
}


template<typename pp>
long long profile_pairing_measurement(
    test_instances_t<pp> elements)
{
    long long start_time = get_nsec_time();

    std::vector<GT<pp>> answers;
    for (size_t i = 0; i < elements.size(); i++) {
        answers.push_back(pp::reduced_pairing(elements[i].first, elements[i].second));
    }

    long long time_delta = get_nsec_time() - start_time;

    return time_delta;
}


template<typename pp>
void print_performance_csv(
    size_t n)
{
    long long time_measurement;
    test_instances_t<pp> elements;

    elements = generate_group_elements<pp>(n);
    time_measurement = profile_pairing_measurement<pp>(elements);

    printf("\t%lld", time_measurement); fflush(stdout);
    printf("\n");



}




template<typename pp>
void profile_pairing(){
    pp::init_public_params();
    print_performance_csv<pp>(100);
}

int main()
{
    inhibit_profiling_info = true;
    inhibit_profiling_counters = true;

    print_compilation_info();

    printf("Profiling BN254\n");
    profile_pairing<bn254_pp>();

    printf("Profiling BN183\n");
    profile_pairing<bn183_pp>();

    printf("Profiling BN124\n");
    profile_pairing<bn124_pp>();

    printf("Profiling EDWARDS181\n");
    profile_pairing<edwards181_pp>();

    printf("Profiling EDWARDS97\n");
    profile_pairing<edwards97_pp>();

    printf("Profiling EDWARDS58\n");
    profile_pairing<edwards58_pp>();

    return 0;
}
