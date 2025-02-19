/**
 *****************************************************************************
 Profile pairings.

 *****************************************************************************/
#include <cstdio>
#include <vector>

#include "libff/algebra/curves/bn254/bn254_pp.hpp"
#include "libff/algebra/curves/bn183/bn183_pp.hpp"
#include "libff/algebra/curves/bn124/bn124_pp.hpp"
#include "libff/algebra/curves/edwards181/edwards181_pp.hpp"
#include "libff/algebra/curves/edwards97/edwards97_pp.hpp"
#include "libff/algebra/curves/edwards58/edwards58_pp.hpp"
#include "libff/algebra/scalar_multiplication/multiexp.hpp"
#include "libff/common/rng.hpp"

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


struct PairingStepsMeasurements
{
    long long precompute_g1;
    long long precompute_g2;
    long long miller_loop;
    long long final_exponentiation;
};

template<typename pp>
PairingStepsMeasurements profile_pairing_steps_measurement(
    test_instances_t<pp> elements)
{
    PairingStepsMeasurements result;

    long long start_time, end_time;

    std::vector<GT<pp>> miller_result;
    std::vector<GT<pp>> reduced;
    std::vector<G1_precomp<pp>> precomputed_g1;
    std::vector<G2_precomp<pp>> precomputed_g2;

    miller_result.reserve(elements.size());
    reduced.reserve(elements.size());
    precomputed_g1.reserve(elements.size());
    precomputed_g2.reserve(elements.size());

    start_time = get_nsec_time();

    for (size_t i = 0; i < elements.size(); i++)
    {
        precomputed_g1.push_back(pp::precompute_G1(elements[i].first));
    }

    end_time = get_nsec_time();
    result.precompute_g1 = end_time - start_time;
    start_time = get_nsec_time();

    for (size_t i = 0; i < elements.size(); i++)
    {
        precomputed_g2.push_back(pp::precompute_G2(elements[i].second));
    }

    end_time = get_nsec_time();
    result.precompute_g2 = end_time - start_time;
    start_time = get_nsec_time();

    for (size_t i = 0; i < elements.size(); i++)
    {
        miller_result.push_back(pp::miller_loop(precomputed_g1[i], precomputed_g2[i]));
    }

    end_time = get_nsec_time();
    result.miller_loop = end_time - start_time;
    start_time = get_nsec_time();

    for (size_t i = 0; i < elements.size(); i++)
    {
        reduced.push_back(pp::final_exponentiation(miller_result[i]));
    }

    end_time = get_nsec_time();
    result.final_exponentiation = end_time - start_time;
    start_time = get_nsec_time();

    return result;
}


template<typename pp>
void print_performance_csv(
    size_t n)
{
    long long time_measurement;
    PairingStepsMeasurements steps_measurements;
    test_instances_t<pp> elements;

    elements = generate_group_elements<pp>(n);
    time_measurement = profile_pairing_measurement<pp>(elements);

    printf("\tReduced Pairing:\t%lld", time_measurement); fflush(stdout);

    steps_measurements = profile_pairing_steps_measurement<pp>(elements);
    printf("\tPrecompute G1:\t%lld", steps_measurements.precompute_g1);
    printf("\tPrecompute G2:\t%lld", steps_measurements.precompute_g2);
    printf("\tMiller Loop:\t%lld", steps_measurements.miller_loop);
    printf("\tFinal Exponentiation:\t%lld", steps_measurements.final_exponentiation);

    fflush(stdout);

    printf("\n");

}




template<typename pp>
void profile_pairing(){
    pp::init_public_params();
    print_performance_csv<pp>(1000);
}

int main()
{
    inhibit_profiling_info = true;
    inhibit_profiling_counters = true;

    print_compilation_info();

    printf("Profiling BN254\n");
    profile_pairing<bn254_pp>();

    /*printf("Profiling BN183\n");
    profile_pairing<bn183_pp>();

    printf("Profiling BN124\n");
    profile_pairing<bn124_pp>();

    printf("Profiling EDWARDS181\n");
    profile_pairing<edwards181_pp>();

    printf("Profiling EDWARDS97\n");
    profile_pairing<edwards97_pp>();*/

    printf("Profiling EDWARDS58\n");
    profile_pairing<edwards58_pp>();

    printf("Profiling EDWARDS181\n");
    profile_pairing<edwards181_pp>();

    return 0;
}
