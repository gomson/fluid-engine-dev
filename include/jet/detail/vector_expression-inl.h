// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_DETAIL_VECTOR_EXPRESSION_INL_H
#define INCLUDE_JET_DETAIL_VECTOR_EXPRESSION_INL_H

#include <jet/vector_expression.h>

namespace jet {

// MARK: VectorExpression

template <typename T, typename E>
size_t VectorExpression<T, E>::size() const {
    return static_cast<const E &>(*this).size();
}

template <typename T, typename E>
const E &VectorExpression<T, E>::operator()() const {
    return static_cast<const E &>(*this);
}

// MARK: VectorUnaryOp

template <typename T, typename E, typename Op>
VectorUnaryOp<T, E, Op>::VectorUnaryOp(const E &u) : _u(u) {}

template <typename T, typename E, typename Op>
size_t VectorUnaryOp<T, E, Op>::size() const {
    return _u.size();
}

template <typename T, typename E, typename Op>
T VectorUnaryOp<T, E, Op>::operator[](size_t i) const {
    return _op(_u[i]);
}

// MARK: VectorBinaryOp

template <typename T, typename E1, typename E2, typename Op>
VectorBinaryOp<T, E1, E2, Op>::VectorBinaryOp(const E1 &u, const E2 &v)
    : _u(u), _v(v) {
    JET_ASSERT(u.size() == v.size());
}

template <typename T, typename E1, typename E2, typename Op>
size_t VectorBinaryOp<T, E1, E2, Op>::size() const {
    return _v.size();
}

template <typename T, typename E1, typename E2, typename Op>
T VectorBinaryOp<T, E1, E2, Op>::operator[](size_t i) const {
    return _op(_u[i], _v[i]);
}

template <typename T, typename E, typename Op>
VectorScalarBinaryOp<T, E, Op>::VectorScalarBinaryOp(const E &u, const T &v)
    : _u(u), _v(v) {}

template <typename T, typename E, typename Op>
size_t VectorScalarBinaryOp<T, E, Op>::size() const {
    return _u.size();
}

template <typename T, typename E, typename Op>
T VectorScalarBinaryOp<T, E, Op>::operator[](size_t i) const {
    return _op(_u[i], _v);
}

}  // namespace jet

#endif  // INCLUDE_JET_DETAIL_VECTOR_EXPRESSION_INL_H
