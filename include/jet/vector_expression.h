// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_VECTOR_EXPRESSION_H_
#define INCLUDE_JET_VECTOR_EXPRESSION_H_

namespace jet {

// MARK: VectorExpression

template <typename T, typename E>
class VectorExpression {
 public:
    size_t size() const;

    const E& operator()() const;
};

// MARK: VectorUnaryOp

template <typename T, typename E, typename Op>
class VectorUnaryOp : public VectorExpression<T, VectorUnaryOp<T, E, Op>> {
 public:
    VectorUnaryOp(const E& u);

    size_t size() const;

    T operator[](size_t i) const;

 private:
    const E& _u;
    Op _op;
};

// MARK: VectorUnaryOp Aliases

template <typename T, typename E, typename U>
using VectorTypeCast = VectorUnaryOp<T, E, TypeCast<U, T>>;

// MARK: VectorBinaryOp

template <typename T, typename E1, typename E2, typename Op>
class VectorBinaryOp
    : public VectorExpression<T, VectorBinaryOp<T, E1, E2, Op>> {
 public:
    VectorBinaryOp(const E1& u, const E2& v);

    size_t size() const;

    T operator[](size_t i) const;

 private:
    const E1& _u;
    const E2& _v;
    Op _op;
};

template <typename T, typename E, typename Op>
class VectorScalarBinaryOp
    : public VectorExpression<T, VectorScalarBinaryOp<T, E, Op>> {
 public:
    VectorScalarBinaryOp(const E& u, const T& v);

    size_t size() const;

    T operator[](size_t i) const;

 private:
    const E& _u;
    T _v;
    Op _op;
};

// MARK: VectorBinaryOp Aliases

template <typename T, typename E1, typename E2>
using VectorAdd = VectorBinaryOp<T, E1, E2, std::plus<T>>;

template <typename T, typename E>
using VectorScalarAdd = VectorScalarBinaryOp<T, E, std::plus<T>>;

template <typename T, typename E1, typename E2>
using VectorSub = VectorBinaryOp<T, E1, E2, std::minus<T>>;

template <typename T, typename E>
using VectorScalarSub = VectorScalarBinaryOp<T, E, std::minus<T>>;

template <typename T, typename E>
using VectorScalarRSub = VectorScalarBinaryOp<T, E, RMinus<T>>;

template <typename T, typename E1, typename E2>
using VectorMul = VectorBinaryOp<T, E1, E2, std::multiplies<T>>;

template <typename T, typename E>
using VectorScalarMul = VectorScalarBinaryOp<T, E, std::multiplies<T>>;

template <typename T, typename E1, typename E2>
using VectorDiv = VectorBinaryOp<T, E1, E2, std::divides<T>>;

template <typename T, typename E>
using VectorScalarDiv = VectorScalarBinaryOp<T, E, std::divides<T>>;

template <typename T, typename E>
using VectorScalarRDiv = VectorScalarBinaryOp<T, E, RDivides<T>>;

// MARK: Global Functions

template <typename T, typename E>
VectorScalarAdd<T, E> operator+(const T& a, const VectorExpression<T, E>& b);

template <typename T, typename E>
VectorScalarAdd<T, E> operator+(const VectorExpression<T, E>& a, const T& b);

template <typename T, typename E1, typename E2>
VectorAdd<T, E1, E2> operator+(const VectorExpression<T, E1>& a,
                               const VectorExpression<T, E2>& b);

template <typename T, typename E>
VectorScalarRSub<T, E> operator-(const T& a, const VectorExpression<T, E>& b);

template <typename T, typename E>
VectorScalarSub<T, E> operator-(const VectorExpression<T, E>& a, const T& b);

template <typename T, typename E1, typename E2>
VectorSub<T, E1, E2> operator-(const VectorExpression<T, E1>& a,
                               const VectorExpression<T, E2>& b);

template <typename T, typename E>
VectorScalarMul<T, E> operator*(const T& a, const VectorExpression<T, E>& b);

template <typename T, typename E>
VectorScalarMul<T, E> operator*(const VectorExpression<T, E>& a, const T& b);

template <typename T, typename E1, typename E2>
VectorMul<T, E1, E2> operator*(const VectorExpression<T, E1>& a,
                               const VectorExpression<T, E2>& b);

template <typename T, typename E>
VectorScalarRDiv<T, E> operator/(const T& a, const VectorExpression<T, E>& b);

template <typename T, typename E>
VectorScalarDiv<T, E> operator/(const VectorExpression<T, E>& a, const T& b);

template <typename T, typename E1, typename E2>
VectorDiv<T, E1, E2> operator/(const VectorExpression<T, E1>& a,
                               const VectorExpression<T, E2>& b);

}  // namespace jet

#include "detail/vector_expression-inl.h"

#endif  // INCLUDE_JET_VECTOR_EXPRESSION_H_
