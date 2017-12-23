#pragma once

#include <array>
#include <algorithm>

template <size_t NUM_DERIVATIVES>
class CDualNumber
{
public:
    CDualNumber(float value = 0.0f)
    {
        m_value = value;
        std::fill(m_derivatives.begin(), m_derivatives.end(), 0.0f);
    }

    float m_value;
    std::array<float, NUM_DERIVATIVES> m_derivatives;
};

// ret.real = a.real + b.real
// ret.dual = a.dual + b.dual
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> operator + (const CDualNumber<NUM_DERIVATIVES> &a, const CDualNumber<NUM_DERIVATIVES> &b)
{
    CDualNumber<NUM_DERIVATIVES> ret;
    ret.m_value = a.m_value + b.m_value;
    for (size_t i = 0; i < NUM_DERIVATIVES; ++i)
        ret.m_derivatives[i] = a.m_derivatives[i] + b.m_derivatives[i];
    return ret;
}

// ret.real = a.real - b.real
// ret.dual = a.dual - b.dual
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> operator - (const CDualNumber<NUM_DERIVATIVES> &a, const CDualNumber<NUM_DERIVATIVES> &b)
{
    CDualNumber<NUM_DERIVATIVES> ret;
    ret.m_value = a.m_value - b.m_value;
    for (size_t i = 0; i < NUM_DERIVATIVES; ++i)
        ret.m_derivatives[i] = a.m_derivatives[i] - b.m_derivatives[i];
    return ret;
}

// ret.real = a.real * b.real
// ret.dual = a.real * b.dual + a.dual * b.real
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> operator * (const CDualNumber<NUM_DERIVATIVES> &a, const CDualNumber<NUM_DERIVATIVES> &b)
{
    CDualNumber<NUM_DERIVATIVES> ret;

    ret.m_value = a.m_value * b.m_value;

    for (size_t i = 0; i < NUM_DERIVATIVES; ++i)
        ret.m_derivatives[i] = a.m_value * b.m_derivatives[i] + a.m_derivatives[i] * b.m_value;

    return ret;
}

// ret.real = a.real / b.real
// ret.dual = (a.dual * b.real - a.real * b.dual) / (b.real * b.real)
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> operator / (const CDualNumber<NUM_DERIVATIVES> &a, const CDualNumber<NUM_DERIVATIVES> &b)
{
    CDualNumber<NUM_DERIVATIVES> ret;

    ret.m_value = a.m_value / b.m_value;

    float bRealSquared = b.m_value * b.m_value;

    for (size_t i = 0; i < NUM_DERIVATIVES; ++i)
        ret.m_derivatives[i] = (a.m_derivatives[i] * b.m_value - a.m_value * b.m_derivatives[i]) / bRealSquared;

    return ret;
}

// ret.real = atan(a.real)
// ret.dual = a.dual / (1 + a.real * a.real)
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> atan (const CDualNumber<NUM_DERIVATIVES> &a)
{
    CDualNumber<NUM_DERIVATIVES> ret;
    ret.m_value = std::atanf(a.m_value);

    float aRealSquared = a.m_value * a.m_value;

    for (size_t i = 0; i < NUM_DERIVATIVES; ++i)
    {
        ret.m_derivatives[i] = a.m_derivatives[i] / (1.0f + aRealSquared);
    }

    return ret;
}

// Lerp done algebraicly
template <size_t NUM_DERIVATIVES>
inline CDualNumber<NUM_DERIVATIVES> Lerp (const CDualNumber<NUM_DERIVATIVES> &a, const CDualNumber<NUM_DERIVATIVES> &b, const CDualNumber<NUM_DERIVATIVES> &t)
{
    return a * (CDualNumber<NUM_DERIVATIVES>(1.0f) - t) + b * t;
}