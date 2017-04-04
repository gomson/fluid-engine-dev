// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_MATRIX_EXPRESSION_H_
#define INCLUDE_JET_MATRIX_EXPRESSION_H_

#include <jet/vector_expression.h>

namespace jet {

// MARK: MatrixExpression

template <typename T, typename E>
class MatrixExpression {
 public:
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    const E& operator()() const;
};

template <typename T>
class MatrixConstant : public MatrixExpression<T, MatrixConstant<T>> {
 public:
    MatrixConstant(size_t m, size_t n, const T& c);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    size_t _m;
    size_t _n;
    T _c;
};

template <typename T>
class MatrixIdentity : public MatrixExpression<T, MatrixIdentity<T>> {
 public:
    MatrixIdentity(size_t m);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    size_t _m;
};

// MARK: MatrixUnaryOp

template <typename T, typename E, typename Op>
class MatrixUnaryOp : public MatrixExpression<T, MatrixUnaryOp<T, E, Op>> {
 public:
    MatrixUnaryOp(const E& u);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E& _u;
    Op _op;
};

template <typename T, typename E>
class MatrixDiagonal : public MatrixExpression<T, MatrixDiagonal<T, E>> {
 public:
    MatrixDiagonal(const E& u, bool isDiag);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E& _u;
    bool _isDiag;
};

template <typename T, typename E>
class MatrixTriangular : public MatrixExpression<T, MatrixTriangular<T, E>> {
 public:
    MatrixTriangular(const E& u, bool isUpper, bool isStrict);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E& _u;
    bool _isUpper;
    bool _isStrict;
};

// MARK: MatrixUnaryOp Aliases

template <typename T, typename E, typename U>
using MatrixTypeCast = MatrixUnaryOp<T, E, TypeCast<U, T>>;

// MARK: MatrixBinaryOp

template <typename T, typename E1, typename E2, typename Op>
class MatrixBinaryOp
    : public MatrixExpression<T, MatrixBinaryOp<T, E1, E2, Op>> {
 public:
    MatrixBinaryOp(const E1& u, const E2& v);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E1& _u;
    const E2& _v;
    Op _op;
};

template <typename T, typename E, typename Op>
class MatrixScalarBinaryOp
    : public MatrixExpression<T, MatrixScalarBinaryOp<T, E, Op>> {
 public:
    MatrixScalarBinaryOp(const E& u, const T& v);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E& _u;
    T _v;
    Op _op;
};

//!
//! \tparam T   Element value type.
//! \tparam ME  Matrix expression.
//! \tparam VE  Vector expression.
//!
template <typename T, typename ME, typename VE>
class MatrixVectorMul : public VectorExpression<T, MatrixVectorMul<T, ME, VE>> {
 public:
    MatrixVectorMul(const ME& m, const VE& v);
    size_t size() const;
    T operator[](size_t i) const;

 private:
    const ME& _m;
    const VE& _v;
};

template <typename T, typename E1, typename E2>
class MatrixMul : public MatrixExpression<T, MatrixMul<T, E1, E2>> {
 public:
    MatrixMul(const E1& u, const E2& v);
    Size2 size() const;
    size_t rows() const;
    size_t cols() const;
    T operator()(size_t i, size_t j) const;

 private:
    const E1& _u;
    const E2& _v;
};

// MARK: MatrixBinaryOp Aliases

template <typename T, typename E1, typename E2>
using MatrixAdd = MatrixBinaryOp<T, E1, E2, std::plus<T>>;

template <typename T, typename E>
using MatrixScalarAdd = MatrixScalarBinaryOp<T, E, std::plus<T>>;

template <typename T, typename E1, typename E2>
using MatrixSub = MatrixBinaryOp<T, E1, E2, std::minus<T>>;

template <typename T, typename E>
using MatrixScalarSub = MatrixScalarBinaryOp<T, E, std::minus<T>>;

template <typename T, typename E>
using MatrixScalarRSub = MatrixScalarBinaryOp<T, E, RMinus<T>>;

template <typename T, typename E>
using MatrixScalarMul = MatrixScalarBinaryOp<T, E, std::multiplies<T>>;

template <typename T, typename E>
using MatrixScalarDiv = MatrixScalarBinaryOp<T, E, std::divides<T>>;

template <typename T, typename E>
using MatrixScalarRDiv = MatrixScalarBinaryOp<T, E, RDivides<T>>;

}  // namespace jet

#include "detail/matrix_expression-inl.h"

#endif  // INCLUDE_JET_MATRIX_EXPRESSION_H_
