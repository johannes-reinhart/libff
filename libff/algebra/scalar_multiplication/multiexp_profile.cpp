/**
 *****************************************************************************
 Implementation of interfaces for multi-exponentiation routines.

 See multiexp.hpp .
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <cstdio>
#include <vector>

#include <libff/algebra/curves/bn254/bn254_pp.hpp>
#include <libff/algebra/curves/bn183/bn183_pp.hpp>
#include <libff/algebra/curves/bn124/bn124_pp.hpp>
#include <libff/algebra/curves/edwards181/edwards181_pp.hpp>
#include <libff/algebra/curves/edwards97/edwards97_pp.hpp>
#include <libff/algebra/curves/edwards58/edwards58_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/rng.hpp>

using namespace libff;

using std::size_t;

template <typename GroupT>
using run_result_t = std::pair<long long, std::vector<GroupT> >;

template <typename T>
using test_instances_t = std::vector<std::vector<T> >;

template<typename GroupT>
test_instances_t<GroupT> generate_group_elements(size_t count, size_t size)
{
    // generating a random group element is expensive,
    // so for now we only generate a single one and repeat it
    test_instances_t<GroupT> result(count);

    for (size_t i = 0; i < count; i++) {
        GroupT x = GroupT::random_element();
#ifdef USE_MIXED_ADDITION
        x.to_special(); // djb requires input to be in special form
#endif
        for (size_t j = 0; j < size; j++) {
            result[i].push_back(x);
            // result[i].push_back(GroupT::random_element());
        }
    }

    return result;
}

template<typename FieldT>
test_instances_t<FieldT> generate_scalars(size_t count, size_t size)
{
    // we use SHA512_rng because it is much faster than
    // FieldT::random_element()
    test_instances_t<FieldT> result(count);

    for (size_t i = 0; i < count; i++) {
        for (size_t j = 0; j < size; j++) {
            result[i].push_back(SHA512_rng<FieldT>(i * size + j));
        }
    }

    return result;
}

template<typename GroupT, typename FieldT, multi_exp_method Method>
run_result_t<GroupT> profile_multiexp(
    test_instances_t<GroupT> group_elements,
    test_instances_t<FieldT> scalars)
{
    long long start_time = get_nsec_time();

    std::vector<GroupT> answers;
    for (size_t i = 0; i < group_elements.size(); i++) {
        answers.push_back(multi_exp<GroupT, FieldT, Method>(
            group_elements[i].cbegin(), group_elements[i].cend(),
            scalars[i].cbegin(), scalars[i].cend(),
            1));
    }

    long long time_delta = get_nsec_time() - start_time;

    return run_result_t<GroupT>(time_delta, answers);
}


template<typename GroupT, typename FieldT>
void print_performance_csv(
    size_t expn_start,
    size_t expn_end_fast,
    size_t expn_end_naive,
    bool compare_answers)
{
    for (size_t expn = expn_start; expn <= expn_end_fast; expn++) {
        printf("%ld", expn); fflush(stdout);

        test_instances_t<GroupT> group_elements =
            generate_group_elements<GroupT>(10, 1 << expn);
        test_instances_t<FieldT> scalars =
            generate_scalars<FieldT>(10, 1 << expn);

        run_result_t<GroupT> result_bos_coster =
            profile_multiexp<GroupT, FieldT, multi_exp_method_bos_coster>(
                group_elements, scalars);
        printf("\t%lld", result_bos_coster.first); fflush(stdout);

        run_result_t<GroupT> result_djb =
            profile_multiexp<GroupT, FieldT, multi_exp_method_BDLO12>(
                group_elements, scalars);
        printf("\t%lld", result_djb.first); fflush(stdout);

        if (compare_answers && (result_bos_coster.second != result_djb.second)) {
            fprintf(stderr, "Answers NOT MATCHING (bos coster != djb)\n");
        }

        if (expn <= expn_end_naive) {
            run_result_t<GroupT> result_naive =
                profile_multiexp<GroupT, FieldT, multi_exp_method_naive>(
                    group_elements, scalars);
            printf("\t%lld", result_naive.first); fflush(stdout);

            if (compare_answers && (result_bos_coster.second != result_naive.second)) {
                fprintf(stderr, "Answers NOT MATCHING (bos coster != naive)\n");
            }
        }

        printf("\n");
    }
}


template<typename GroupT, typename FieldT>
void field_print_performance_csv(
        size_t expn_start,
        size_t expn_end_fast,
        size_t expn_end_naive,
        bool compare_answers)
{
    for (size_t expn = expn_start; expn <= expn_end_fast; expn++) {
        printf("%ld", expn); fflush(stdout);

        test_instances_t<GroupT> group_elements =
                generate_group_elements<GroupT>(10, 1 << expn);
        test_instances_t<FieldT> scalars =
                generate_scalars<FieldT>(10, 1 << expn);

        if (expn <= expn_end_naive) {
            run_result_t<GroupT> result_naive =
                    profile_multiexp<GroupT, FieldT, multi_exp_method_naive_plain>(
                            group_elements, scalars);
            printf("\t%lld", result_naive.first); fflush(stdout);

        }

        printf("\n");
    }
}

template<typename pp>
void profile_groups(){
    printf("G1\n");
    pp::init_public_params();
    print_performance_csv<G1<pp>, Fr<pp> >(2, 12, 12, true);

    printf("G2\n");
    print_performance_csv<G2<pp>, Fr<pp> >(2, 12, 12, true);

    printf("GT\n");
    field_print_performance_csv<GT<pp>, Fq<pp> >(2, 12, 12, true);
}

int main()
{
    print_compilation_info();

    printf("Profiling BN254\n");
    profile_groups<bn254_pp>();

    printf("Profiling BN183\n");
    profile_groups<bn183_pp>();

    printf("Profiling BN124\n");
    profile_groups<bn124_pp>();

    printf("Profiling EDWARDS181\n");
    profile_groups<edwards181_pp>();

    printf("Profiling EDWARDS97\n");
    profile_groups<edwards97_pp>();

    printf("Profiling EDWARDS58\n");
    profile_groups<edwards58_pp>();

    return 0;
}
