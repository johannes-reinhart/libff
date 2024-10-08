/** @file
 *****************************************************************************
 Implementation of misc. math and serialization utility functions
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FIELD_UTILS_TCC_
#define FIELD_UTILS_TCC_

#include <complex>
#include <stdexcept>

namespace libff {

using std::size_t;

template<typename FieldT>
field_type get_field_type(const typename enable_if<is_multiplicative<FieldT>::value, FieldT>::type elem)
{
    UNUSED(elem); // only to identify field type
    return multiplicative_field_type;
}

template<typename FieldT>
field_type get_field_type(const typename enable_if<is_additive<FieldT>::value, FieldT>::type elem)
{
    UNUSED(elem); // only to identify field type
    return additive_field_type;
}

template<typename FieldT>
std::size_t log_of_field_size_helper(
    typename enable_if<is_multiplicative<FieldT>::value, FieldT>::type field_elem)
{
    UNUSED(field_elem);
    return FieldT::ceil_size_in_bits();
}

template<typename FieldT>
std::size_t log_of_field_size_helper(
    typename enable_if<is_additive<FieldT>::value, FieldT>::type field_elem)
{
    UNUSED(field_elem);
    return FieldT::extension_degree();
}

template<typename FieldT>
std::size_t soundness_log_of_field_size_helper(
    typename enable_if<is_multiplicative<FieldT>::value, FieldT>::type field_elem)
{
    UNUSED(field_elem);
    /** size in bits is the number of bits needed to represent a field element.
     *  However there isn't perfect alignment between the number of bits and the number of field elements,
     *  there could be a factor of two difference.
     *  For calculating soundness, we use the log of field size as number of bits - 1,
     *  as (2 << returned) size lower bounds the actual size.
    */
    return FieldT::ceil_size_in_bits() - 1;
}

template<typename FieldT>
std::size_t soundness_log_of_field_size_helper(
    typename enable_if<is_additive<FieldT>::value, FieldT>::type field_elem)
{
    UNUSED(field_elem);
    return FieldT::extension_degree();
}

template<typename FieldT>
std::size_t get_word_of_field_elem(
    typename enable_if<is_additive<FieldT>::value, FieldT>::type field_elem, size_t word)
{
    return field_elem.to_words()[word];
}

template<typename FieldT>
std::size_t get_word_of_field_elem(
    typename enable_if<is_multiplicative<FieldT>::value, FieldT>::type field_elem, size_t word)
{
    return field_elem.as_bigint().data[word];
}

template<typename FieldT>
FieldT coset_shift()
{
    return FieldT::multiplicative_generator.squared();
}

template<typename FieldT>
typename std::enable_if<std::is_same<FieldT, Double>::value, FieldT>::type
get_root_of_unity(const size_t n)
{
    return FieldT(cos(2 * PI / n), sin(2 * PI / n));
}

template<typename FieldT>
typename std::enable_if<!std::is_same<FieldT, Double>::value, FieldT>::type
get_root_of_unity(const size_t n)
{
    const size_t logn = log2(n);
    if (n != (1u << logn)) throw std::invalid_argument("libff::get_root_of_unity: expected n == (1u << logn)");
    if (logn > FieldT::s) throw std::invalid_argument("libff::get_root_of_unity: expected logn <= FieldT::s");

    FieldT omega = FieldT::root_of_unity;
    for (size_t i = FieldT::s; i > logn; --i)
    {
        omega *= omega;
    }

    return omega;
}

template<typename FieldT>
std::vector<FieldT> pack_int_vector_into_field_element_vector(const std::vector<size_t> &v, const size_t w)
{
    const size_t chunk_bits = FieldT::floor_size_in_bits();
    const size_t repacked_size = div_ceil(v.size() * w, chunk_bits);
    std::vector<FieldT> result(repacked_size);

    for (size_t i = 0; i < repacked_size; ++i)
    {
        bigint<FieldT::num_limbs> b;
        for (size_t j = 0; j < chunk_bits; ++j)
        {
            const size_t word_index = (i * chunk_bits + j) / w;
            const size_t pos_in_word = (i * chunk_bits + j) % w;
            const size_t word_or_0 = (word_index < v.size() ? v[word_index] : 0);
            const size_t bit = (word_or_0 >> pos_in_word) & 1;

            b.data[j / GMP_NUMB_BITS] |= bit << (j % GMP_NUMB_BITS);
        }
        result[i] = FieldT(b);
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> pack_bit_vector_into_field_element_vector(const bit_vector &v, const size_t chunk_bits)
{
    assert(chunk_bits <= FieldT::floor_size_in_bits());

    const size_t repacked_size = div_ceil(v.size(), chunk_bits);
    std::vector<FieldT> result(repacked_size);

    for (size_t i = 0; i < repacked_size; ++i)
    {
        bigint<FieldT::num_limbs> b;
        for (size_t j = 0; j < chunk_bits; ++j)
        {
            b.data[j / GMP_NUMB_BITS] |= ((i * chunk_bits + j) < v.size() && v[i * chunk_bits + j] ? 1ll : 0ll) << (j % GMP_NUMB_BITS);
        }
        result[i] = FieldT(b);
    }

    return result;
}

template<typename FieldT>
std::vector<FieldT> pack_bit_vector_into_field_element_vector(const bit_vector &v)
{
    return pack_bit_vector_into_field_element_vector<FieldT>(v, FieldT::floor_size_in_bits());
}

template<typename FieldT>
std::vector<FieldT> convert_bit_vector_to_field_element_vector(const bit_vector &v)
{
    std::vector<FieldT> result;
    result.reserve(v.size());

    for (const bool b : v)
    {
        result.emplace_back(b ? FieldT::one() : FieldT::zero());
    }

    return result;
}

template<typename FieldT>
bit_vector convert_field_element_vector_to_bit_vector(const std::vector<FieldT> &v)
{
    bit_vector result;

    for (const FieldT &el : v)
    {
        const bit_vector el_bits = convert_field_element_to_bit_vector<FieldT>(el);
        result.insert(result.end(), el_bits.begin(), el_bits.end());
    }

    return result;
}

template<typename FieldT>
bit_vector convert_field_element_to_bit_vector(const FieldT &el)
{
    bit_vector result;

    bigint<FieldT::num_limbs> b = el.as_bigint();
    for (size_t i = 0; i < FieldT::ceil_size_in_bits(); ++i)
    {
        result.push_back(b.test_bit(i));
    }

    return result;
}

template<typename FieldT>
bit_vector convert_field_element_to_bit_vector(const FieldT &el, const size_t bitcount)
{
    bit_vector result = convert_field_element_to_bit_vector(el);
    result.resize(bitcount);

    return result;
}

template<typename FieldT>
FieldT convert_bit_vector_to_field_element(const bit_vector &v)
{
    assert(v.size() <= FieldT::ceil_size_in_bits());

    FieldT res = FieldT::zero();
    FieldT c = FieldT::one();
    for (bool b : v)
    {
        res += b ? c : FieldT::zero();
        c += c;
    }
    return res;
}

template<typename FieldT>
void batch_invert(std::vector<FieldT> &vec)
{
    std::vector<FieldT> prod;
    prod.reserve(vec.size());

    FieldT acc = FieldT::one();

    for (auto el : vec)
    {
        assert(!el.is_zero());
        prod.emplace_back(acc);
        acc = acc * el;
    }

    FieldT acc_inverse = acc.inverse();

    for (long i = static_cast<long>(vec.size()-1); i >= 0; --i)
    {
        const FieldT old_el = vec[i];
        vec[i] = acc_inverse * prod[i];
        acc_inverse = acc_inverse * old_el;
    }
}

template<typename FieldT>
bool is_negative(FieldT v){
    // ethsnarks utils has it the other way around
    // According to Original Bernstein EdDSA Paper:
    // "We use the encoding of Fq to define some field elements as being negative: specifically,
    // x is negative if the (b−1)-bit encoding of x is lexicographically larger than the (b−1)-bit encoding of −x"
    // -> Depends on endianness (?)
    // -> Here we compare numbers -> lexicographically big-endian
    FieldT mv = -v;
    return mv.as_bigint() < v.as_bigint();
}

} // namespace libff
#endif // FIELD_UTILS_TCC_
