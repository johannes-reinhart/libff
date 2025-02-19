/** @file
 *****************************************************************************
 Implementation of serialization routines.

 See serialization.hpp .
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef SERIALIZATION_TCC_
#define SERIALIZATION_TCC_

#include <cassert>
#include <sstream>
#include <iomanip>

#include <libff/common/utils.hpp>

namespace libff {

using std::size_t;

inline void consume_newline(std::istream &in)
{
    char c;
    in.read(&c, 1);
}

inline void consume_OUTPUT_NEWLINE(std::istream &in)
{
#ifdef BINARY_OUTPUT
    // nothing to consume
    UNUSED(in);
#else
    char c;
    in.read(&c, 1);
#endif
}

inline void consume_OUTPUT_SEPARATOR(std::istream &in)
{
#ifdef BINARY_OUTPUT
    // nothing to consume
    UNUSED(in);
#else
    char c;
    in.read(&c, 1);
#endif
}

inline void output_bool(std::ostream &out, const bool b)
{
    out << (b ? 1 : 0) << "\n";
}

inline void input_bool(std::istream &in, bool &b)
{
    size_t tmp;
    in >> tmp;
    consume_newline(in);
    assert(tmp == 0 || tmp == 1);

    b = (tmp == 1 ? true : false);
}

inline void output_bool_vector(std::ostream &out, const std::vector<bool> &v)
{
    out << v.size() << "\n";
    for (const bool b : v)
    {
        output_bool(out, b);
    }
}

inline void input_bool_vector(std::istream &in, std::vector<bool> &v)
{
    size_t size;
    in >> size;
    consume_newline(in);
    v.resize(size);
    for (size_t i = 0; i < size; ++i)
    {
        bool b;
        input_bool(in, b);
        v[i] = b;
    }
}

inline void output_bytes(std::ostream& out, const std::vector<uint8_t> &v)
{
#ifdef BINARY_OUTPUT
    assert(v.size() < 256);
    auto size = static_cast<uint8_t>(v.size());
    out.write(reinterpret_cast<const char*>(&size), 1);
    out.write(reinterpret_cast<const char*>(v.data()), size);
#else
    out << std::hex << std::setfill('0');
    for (const auto& byte : v) {
        out << std::setw(2) << static_cast<int>(byte);
    }
    out << std::dec << std::setfill(' ');  // Reset to defaults after operation
    out << OUTPUT_SEPARATOR;
#endif
}

inline void input_bytes(std::istream& in, std::vector<uint8_t> &v)
{
#ifdef BINARY_OUTPUT
    uint8_t size;
    in.read(reinterpret_cast<char*>(&size), 1);
    v.clear();
    v.insert(v.end(), size, 0);
    in.read(reinterpret_cast<char*>(v.data()), size);
#else
    std::string s;
    in >> s;
    libff::consume_OUTPUT_SEPARATOR(in);

    assert(s.length() % 2 == 0);
    v.reserve(s.length() / 2);

    // Convert each pair of hexadecimal digits to a byte
    for (size_t i = 0; i < s.length(); i += 2) {
        // Extract two characters from the string
        char highNibble = s[i];
        char lowNibble = s[i + 1];

        // Ensure they are valid hexadecimal digits
        assert(std::isxdigit(highNibble) && std::isxdigit(lowNibble));

        // Convert hex pair to a byte
        uint8_t byte = (std::stoi(s.substr(i, 2), nullptr, 16) & 0xFF);
        v.push_back(byte);
    }
#endif
}

template<typename T>
T reserialize(const T &obj)
{
    std::stringstream ss;
    ss << obj;
    T tmp;
    ss >> tmp;
    return tmp;
}

template<typename T>
size_t get_serialized_size(const T& obj){
    std::stringstream ss;
    ss << obj;
    return ss.str().size();
}



template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T> &v)
{
    static_assert(!std::is_same<T, bool>::value, "this does not work for std::vector<bool>");
    out << v.size() << "\n";
    for (const T& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

template<typename T>
std::istream& operator>>(std::istream& in, std::vector<T> &v)
{
    static_assert(!std::is_same<T, bool>::value, "this does not work for std::vector<bool>");
    size_t size;
    in >> size;
    consume_newline(in);

    v.resize(0);
    for (size_t i = 0; i < size; ++i)
    {
        T elt;
        in >> elt;
        consume_OUTPUT_NEWLINE(in);
        v.push_back(elt);
    }

    return in;
}

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1, T2> &m)
{
    out << m.size() << "\n";

    for (auto &it : m)
    {
        out << it.first << "\n";
        out << it.second << "\n";
    }

    return out;
}

template<typename T1, typename T2>
std::istream& operator>>(std::istream& in, std::map<T1, T2> &m)
{
    m.clear();
    size_t size;
    in >> size;
    consume_newline(in);

    for (size_t i = 0; i < size; ++i)
    {
        T1 k;
        T2 v;
        in >> k;
        consume_newline(in);
        in >> v;
        consume_newline(in);
        m[k] = v;
    }

    return in;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T> &s)
{
    out << s.size() << "\n";

    for (auto &el : s)
    {
        out << el << "\n";
    }

    return out;
}


template<typename T>
std::istream& operator>>(std::istream& in, std::set<T> &s)
{
    s.clear();
    size_t size;
    in >> size;
    consume_newline(in);

    for (size_t i = 0; i < size; ++i)
    {
        T el;
        in >> el;
        consume_newline(in);
        s.insert(el);
    }

    return in;
}

}

#endif // SERIALIZATION_TCC_
