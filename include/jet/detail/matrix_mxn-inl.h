// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_DETAIL_MATRIX_MXN_INL_H_
#define INCLUDE_JET_DETAIL_MATRIX_MXN_INL_H_

#include <jet/macros.h>
#include <jet/math_utils.h>
#include <jet/matrix_mxn.h>
#include <jet/parallel.h>

namespace jet {

// MARK: MatrixExpression

template <typename T, typename E>
Size2 MatrixExpression<T, E>::size() const {
    return static_cast<const E&>(*this).size();
}

template <typename T, typename E>
size_t MatrixExpression<T, E>::rows() const {
    return static_cast<const E&>(*this).rows();
}

template <typename T, typename E>
size_t MatrixExpression<T, E>::cols() const {
    return static_cast<const E&>(*this).cols();
}

template <typename T, typename E>
const E& MatrixExpression<T, E>::operator()() const {
    return static_cast<const E&>(*this);
}

//

template <typename T>
MatrixConstant<T>::MatrixConstant(size_t m, size_t n, const T& c)
    : _m(m), _n(n), _c(c) {}

template <typename T>
Size2 MatrixConstant<T>::size() const {
    return Size2(rows(), cols());
}

template <typename T>
size_t MatrixConstant<T>::rows() const {
    return _m;
}

template <typename T>
size_t MatrixConstant<T>::cols() const {
    return _n;
}

template <typename T>
T MatrixConstant<T>::operator()(size_t, size_t) const {
    return _c;
}

//

template <typename T>
MatrixIdentity<T>::MatrixIdentity(size_t m) : _m(m) {}

template <typename T>
Size2 MatrixIdentity<T>::size() const {
    return Size2(_m, _m);
}

template <typename T>
size_t MatrixIdentity<T>::rows() const {
    return _m;
}

template <typename T>
size_t MatrixIdentity<T>::cols() const {
    return _m;
}

template <typename T>
T MatrixIdentity<T>::operator()(size_t i, size_t j) const {
    return (i == j) ? 1 : 0;
}

// MARK: MatrixUnaryOp

template <typename T, typename E, typename Op>
MatrixUnaryOp<T, E, Op>::MatrixUnaryOp(const E& u) : _u(u) {}

template <typename T, typename E, typename Op>
Size2 MatrixUnaryOp<T, E, Op>::size() const {
    return _u.size();
}

template <typename T, typename E, typename Op>
size_t MatrixUnaryOp<T, E, Op>::rows() const {
    return _u.rows();
}

template <typename T, typename E, typename Op>
size_t MatrixUnaryOp<T, E, Op>::cols() const {
    return _u.cols();
}

template <typename T, typename E, typename Op>
T MatrixUnaryOp<T, E, Op>::operator()(size_t i, size_t j) const {
    return _op(_u(i, j));
}

//

template <typename T, typename E>
MatrixDiagonal<T, E>::MatrixDiagonal(const E& u, bool isDiag)
    : _u(u), _isDiag(isDiag) {}

template <typename T, typename E>
Size2 MatrixDiagonal<T, E>::size() const {
    return _u.size();
}

template <typename T, typename E>
size_t MatrixDiagonal<T, E>::rows() const {
    return _u.rows();
}

template <typename T, typename E>
size_t MatrixDiagonal<T, E>::cols() const {
    return _u.cols();
}

template <typename T, typename E>
T MatrixDiagonal<T, E>::operator()(size_t i, size_t j) const {
    if (_isDiag) {
        return (i == j) ? _u(i, j) : 0;
    } else {
        return (i != j) ? _u(i, j) : 0;
    }
}

//

template <typename T, typename E>
MatrixTriangular<T, E>::MatrixTriangular(const E& u, bool isUpper,
                                         bool isStrict)
    : _u(u), _isUpper(isUpper), _isStrict(isStrict) {}

template <typename T, typename E>
Size2 MatrixTriangular<T, E>::size() const {
    return _u.size();
}

template <typename T, typename E>
size_t MatrixTriangular<T, E>::rows() const {
    return _u.rows();
}

template <typename T, typename E>
size_t MatrixTriangular<T, E>::cols() const {
    return _u.cols();
}

template <typename T, typename E>
T MatrixTriangular<T, E>::operator()(size_t i, size_t j) const {
    if (i < j) {
        return (_isUpper) ? _u(i, j) : 0;
    } else if (i > j) {
        return (!_isUpper) ? _u(i, j) : 0;
    } else {
        return (!_isStrict) ? _u(i, j) : 0;
    }
}

// MARK: MatrixBinaryOp

template <typename T, typename E1, typename E2, typename Op>
MatrixBinaryOp<T, E1, E2, Op>::MatrixBinaryOp(const E1& u, const E2& v)
    : _u(u), _v(v) {
    JET_ASSERT(u.size() == v.size());
}

template <typename T, typename E1, typename E2, typename Op>
Size2 MatrixBinaryOp<T, E1, E2, Op>::size() const {
    return _v.size();
}

template <typename T, typename E1, typename E2, typename Op>
size_t MatrixBinaryOp<T, E1, E2, Op>::rows() const {
    return _v.rows();
}

template <typename T, typename E1, typename E2, typename Op>
size_t MatrixBinaryOp<T, E1, E2, Op>::cols() const {
    return _v.cols();
}

template <typename T, typename E1, typename E2, typename Op>
T MatrixBinaryOp<T, E1, E2, Op>::operator()(size_t i, size_t j) const {
    return _op(_u(i, j), _v(i, j));
}

// MARK: MatrixScalarBinaryOp

template <typename T, typename E, typename Op>
MatrixScalarBinaryOp<T, E, Op>::MatrixScalarBinaryOp(const E& u, const T& v)
    : _u(u), _v(v) {}

template <typename T, typename E, typename Op>
Size2 MatrixScalarBinaryOp<T, E, Op>::size() const {
    return _u.size();
}

template <typename T, typename E, typename Op>
size_t MatrixScalarBinaryOp<T, E, Op>::rows() const {
    return _u.rows();
}

template <typename T, typename E, typename Op>
size_t MatrixScalarBinaryOp<T, E, Op>::cols() const {
    return _u.cols();
}

template <typename T, typename E, typename Op>
T MatrixScalarBinaryOp<T, E, Op>::operator()(size_t i, size_t j) const {
    return _op(_u(i, j), _v);
}

//

template <typename T, typename ME, typename VE>
MatrixVectorMul<T, ME, VE>::MatrixVectorMul(const ME& m, const VE& v)
    : _m(m), _v(v) {
    JET_ASSERT(_m.cols() == _v.size());
}

template <typename T, typename ME, typename VE>
size_t MatrixVectorMul<T, ME, VE>::size() const {
    return _v.size();
}

template <typename T, typename ME, typename VE>
T MatrixVectorMul<T, ME, VE>::operator[](size_t i) const {
    T sum = 0;
    const size_t n = _m.cols();
    for (size_t j = 0; j < n; ++j) {
        sum += _m(i, j) * _v[j];
    }
    return sum;
}

// MARK: MatrixMul

template <typename T, typename E1, typename E2>
MatrixMul<T, E1, E2>::MatrixMul(const E1& u, const E2& v) : _u(u), _v(v) {
    JET_ASSERT(_u.cols() == _v.rows());
}

template <typename T, typename E1, typename E2>
Size2 MatrixMul<T, E1, E2>::size() const {
    return Size2(_u.rows(), _v.cols());
}

template <typename T, typename E1, typename E2>
size_t MatrixMul<T, E1, E2>::rows() const {
    return _u.rows();
}

template <typename T, typename E1, typename E2>
size_t MatrixMul<T, E1, E2>::cols() const {
    return _v.cols();
}

template <typename T, typename E1, typename E2>
T MatrixMul<T, E1, E2>::operator()(size_t i, size_t j) const {
    // Unoptimized mat-mat-mul
    T sum = 0;
    const size_t n = _u.cols();
    for (size_t k = 0; k < n; ++k) {
        sum += _u(i, k) * _v(k, j);
    }
    return sum;
}

// MARK: MatrixMxN

template <typename T>
MatrixMxN<T>::MatrixMxN() {}

template <typename T>
MatrixMxN<T>::MatrixMxN(size_t m, size_t n, const T& s) {
    resize(m, n, s);
}

template <typename T>
MatrixMxN<T>::MatrixMxN(
    const std::initializer_list<std::initializer_list<T>>& lst) {
    _elements.set(lst);
}

template <typename T>
template <typename E>
MatrixMxN<T>::MatrixMxN(const MatrixExpression<T, E>& other) {
    set(other);
}

template <typename T>
MatrixMxN<T>::MatrixMxN(const MatrixMxN& other) {
    set(other);
}

template <typename T>
MatrixMxN<T>::MatrixMxN(MatrixMxN&& other) {
    (*this) = std::move(other);
}

template <typename T>
MatrixMxN<T>::MatrixMxN(size_t m, size_t n, const T* arr) {
    set(m, n, arr);
}

template <typename T>
void MatrixMxN<T>::resize(size_t m, size_t n, const T& s) {
    // Note that m and n are flipped.
    _elements.resize(n, m, s);
}

template <typename T>
void MatrixMxN<T>::set(const T& s) {
    _elements.set(s);
}

template <typename T>
void MatrixMxN<T>::set(
    const std::initializer_list<std::initializer_list<T>>& lst) {
    _elements.set(lst);
}

template <typename T>
template <typename E>
void MatrixMxN<T>::set(const MatrixExpression<T, E>& other) {
    resize(other.rows(), other.cols());

    // Parallel evaluation of the expression
    const E& expression = other();
    parallelForEachIndex(
        [&](size_t i, size_t j) { (*this)(i, j) = expression(i, j); });
}

template <typename T>
void MatrixMxN<T>::set(size_t m, size_t n, const T* arr) {
    resize(m, n);
    const size_t sz = m * n;
    for (size_t i = 0; i < sz; ++i) {
        _elements[i] = arr[i];
    }
}

template <typename T>
void MatrixMxN<T>::setDiagonal(const T& s) {
    const size_t l = std::min(rows(), cols());
    for (size_t i = 0; i < l; ++i) {
        (*this)(i, i) = s;
    }
}

template <typename T>
void MatrixMxN<T>::setOffDiagonal(const T& s) {
    parallelForEachIndex([&](size_t i, size_t j) {
        if (i != j) {
            (*this)(i, j) = s;
        }
    });
}

template <typename T>
template <typename E>
void MatrixMxN<T>::setRow(size_t i, const VectorExpression<T, E>& row) {
    JET_ASSERT(cols() == row.size());

    const E& e = row();
    parallelFor(kZeroSize, cols(), [&](size_t j) { (*this)(i, j) = e[j]; });
}

template <typename T>
template <typename E>
void MatrixMxN<T>::setColumn(size_t j, const VectorExpression<T, E>& col) {
    JET_ASSERT(rows() == col.size());

    const E& e = col();
    parallelFor(kZeroSize, rows(), [&](size_t i) { (*this)(i, j) = e[i]; });
}

template <typename T>
template <typename E>
bool MatrixMxN<T>::isEqual(const MatrixExpression<T, E>& other) const {
    if (size() != other.size()) {
        return false;
    }

    const E& e = other();
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            if ((*this)(i, j) != e(i, j)) {
                return false;
            }
        }
    }

    return true;
}

template <typename T>
template <typename E>
bool MatrixMxN<T>::isSimilar(const MatrixExpression<T, E>& other,
                             double tol) const {
    if (size() != other.size()) {
        return false;
    }

    const E& e = other();
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            if (std::fabs((*this)(i, j) - e(i, j)) > tol) {
                printf("%f\n", std::fabs((*this)(i, j) - e(i, j)));
                return false;
            }
        }
    }

    return true;
}

template <typename T>
bool MatrixMxN<T>::isSquare() const {
    return rows() == cols();
}

template <typename T>
Size2 MatrixMxN<T>::size() const {
    return Size2(rows(), cols());
}

template <typename T>
size_t MatrixMxN<T>::rows() const {
    return _elements.height();
}

template <typename T>
size_t MatrixMxN<T>::cols() const {
    return _elements.width();
}

template <typename T>
T* MatrixMxN<T>::data() {
    return _elements.data();
}

template <typename T>
const T* const MatrixMxN<T>::data() const {
    return _elements.data();
}

template <typename T>
typename MatrixMxN<T>::Iterator MatrixMxN<T>::begin() {
    return _elements.begin();
}

template <typename T>
typename MatrixMxN<T>::ConstIterator MatrixMxN<T>::begin() const {
    return _elements.begin();
}

template <typename T>
typename MatrixMxN<T>::Iterator MatrixMxN<T>::end() {
    return _elements.end();
}

template <typename T>
typename MatrixMxN<T>::ConstIterator MatrixMxN<T>::end() const {
    return _elements.end();
}

template <typename T>
MatrixScalarAdd<T, MatrixMxN<T>> MatrixMxN<T>::add(const T& s) const {
    return MatrixScalarAdd<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
template <typename E>
MatrixAdd<T, MatrixMxN<T>, E> MatrixMxN<T>::add(const E& m) const {
    return MatrixAdd<T, MatrixMxN, E>(*this, m);
}

template <typename T>
MatrixScalarSub<T, MatrixMxN<T>> MatrixMxN<T>::sub(const T& s) const {
    return MatrixScalarSub<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
template <typename E>
MatrixSub<T, MatrixMxN<T>, E> MatrixMxN<T>::sub(const E& m) const {
    return MatrixSub<T, MatrixMxN, E>(*this, m);
}

template <typename T>
MatrixScalarMul<T, MatrixMxN<T>> MatrixMxN<T>::mul(const T& s) const {
    return MatrixScalarMul<T, MatrixMxN>(*this, s);
}

template <typename T>
template <typename VE>
MatrixVectorMul<T, MatrixMxN<T>, VE> MatrixMxN<T>::mul(
    const VectorExpression<T, VE>& v) const {
    return MatrixVectorMul<T, MatrixMxN<T>, VE>(*this, v());
}

template <typename T>
template <typename E>
MatrixMul<T, MatrixMxN<T>, E> MatrixMxN<T>::mul(const E& m) const {
    return MatrixMul<T, MatrixMxN, E>(*this, m);
}

template <typename T>
MatrixScalarDiv<T, MatrixMxN<T>> MatrixMxN<T>::div(const T& s) const {
    return MatrixScalarDiv<T, MatrixMxN>(*this, s);
}

template <typename T>
MatrixScalarAdd<T, MatrixMxN<T>> MatrixMxN<T>::radd(const T& s) const {
    return MatrixScalarAdd<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
template <typename E>
MatrixAdd<T, MatrixMxN<T>, E> MatrixMxN<T>::radd(const E& m) const {
    return MatrixAdd<T, MatrixMxN<T>, E>(m, *this);
}

template <typename T>
MatrixScalarRSub<T, MatrixMxN<T>> MatrixMxN<T>::rsub(const T& s) const {
    return MatrixScalarRSub<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
template <typename E>
MatrixSub<T, MatrixMxN<T>, E> MatrixMxN<T>::rsub(const E& m) const {
    return MatrixSub<T, MatrixMxN<T>, E>(m, *this);
}

template <typename T>
MatrixScalarMul<T, MatrixMxN<T>> MatrixMxN<T>::rmul(const T& s) const {
    return MatrixScalarMul<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
template <typename E>
MatrixMul<T, MatrixMxN<T>, E> MatrixMxN<T>::rmul(const E& m) const {
    return MatrixMul<T, MatrixMxN<T>, E>(m, *this);
}

template <typename T>
MatrixScalarRDiv<T, MatrixMxN<T>> MatrixMxN<T>::rdiv(const T& s) const {
    return MatrixScalarRDiv<T, MatrixMxN<T>>(*this, s);
}

template <typename T>
void MatrixMxN<T>::iadd(const T& s) {
    set(add(s));
}

template <typename T>
template <typename E>
void MatrixMxN<T>::iadd(const E& m) {
    set(add(m));
}

template <typename T>
void MatrixMxN<T>::isub(const T& s) {
    set(sub(s));
}

template <typename T>
template <typename E>
void MatrixMxN<T>::isub(const E& m) {
    set(sub(m));
}

template <typename T>
void MatrixMxN<T>::imul(const T& s) {
    set(mul(s));
}

template <typename T>
template <typename E>
void MatrixMxN<T>::imul(const E& m) {
    MatrixMxN tmp = mul(m);
    set(tmp);
}

template <typename T>
void MatrixMxN<T>::idiv(const T& s) {
    set(div(s));
}

template <typename T>
void MatrixMxN<T>::transpose() {
    set(transposed());
}

template <typename T>
void MatrixMxN<T>::invert() {
    JET_ASSERT(isSquare());

    // Computes inverse matrix using Gaussian elimination method.
    // https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
    size_t n = rows();
    MatrixMxN& a = *this;
    MatrixMxN rhs = makeIdentity(n);

    for (size_t i = 0; i < n; ++i) {
        // Search for maximum in this column
        T maxEl = std::fabs(a(i, i));
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(a(k, i)) > maxEl) {
                maxEl = std::fabs(a(k, i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        if (maxRow != i) {
            for (size_t k = i; k < n; ++k) {
                std::swap(a(maxRow, k), a(i, k));
                std::swap(rhs(maxRow, k), rhs(i, k));
            }
        }

        // Make all rows except this one 0 in current column
        for (size_t k = 0; k < n; ++k) {
            if (k == i) {
                continue;
            }
            T c = -a(k, i) / a(i, i);
            for (size_t j = 0; j < n; ++j) {
                rhs(k, j) += c * rhs(i, j);
                if (i == j) {
                    a(k, j) = 0;
                } else if (i < j) {
                    a(k, j) += c * a(i, j);
                }
            }
        }

        // Scale this row
        for (size_t i = 0; i < n; ++i) {
            T c = 1 / a(i, i);
            for (size_t j = 0; j < n; ++j) {
                a(i, j) *= c;
                rhs(i, j) *= c;
            }
        }
    }

    set(rhs);
}

template <typename T>
T MatrixMxN<T>::sum() const {
    return parallelReduce(kZeroSize, rows() * cols(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result += _elements[i];
                              }
                              return result;
                          },
                          std::plus<T>());
}

template <typename T>
T MatrixMxN<T>::avg() const {
    return sum() / (rows() * cols());
}

template <typename T>
T MatrixMxN<T>::min() const {
    const T& (*_min)(const T&, const T&) = std::min<T>;
    return parallelReduce(kZeroSize, rows() * cols(),
                          std::numeric_limits<T>::max(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = std::min(result, _elements[i]);
                              }
                              return result;
                          },
                          _min);
}

template <typename T>
T MatrixMxN<T>::max() const {
    const T& (*_max)(const T&, const T&) = std::max<T>;
    return parallelReduce(kZeroSize, rows() * cols(),
                          std::numeric_limits<T>::min(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = std::max(result, _elements[i]);
                              }
                              return result;
                          },
                          _max);
}

template <typename T>
T MatrixMxN<T>::absmin() const {
    return parallelReduce(kZeroSize, rows() * cols(),
                          std::numeric_limits<T>::max(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = jet::absmin(result, _elements[i]);
                              }
                              return result;
                          },
                          jet::absmin<T>);
}

template <typename T>
T MatrixMxN<T>::absmax() const {
    return parallelReduce(kZeroSize, rows() * cols(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = jet::absmax(result, _elements[i]);
                              }
                              return result;
                          },
                          jet::absmax<T>);
}

template <typename T>
T MatrixMxN<T>::trace() const {
    JET_ASSERT(isSquare());
    return parallelReduce(kZeroSize, rows(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result += _elements(i, i);
                              }
                              return result;
                          },
                          std::plus<T>());
}

template <typename T>
T MatrixMxN<T>::determinant() const {
    JET_ASSERT(isSquare());

    // Computes inverse matrix using Gaussian elimination method.
    // https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
    size_t n = rows();
    MatrixMxN a(*this);

    T result = 1;
    for (size_t i = 0; i < n; ++i) {
        // Search for maximum in this column
        T maxEl = std::fabs(a(i, i));
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(a(k, i)) > maxEl) {
                maxEl = std::fabs(a(k, i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        if (maxRow != i) {
            for (size_t k = i; k < n; ++k) {
                std::swap(a(maxRow, k), a(i, k));
            }
            result *= -1;
        }

        // Make all rows below this one 0 in current column
        for (size_t k = i + 1; k < n; ++k) {
            T c = -a(k, i) / a(i, i);
            for (size_t j = i; j < n; ++j) {
                if (i == j) {
                    a(k, j) = 0;
                } else {
                    a(k, j) += c * a(i, j);
                }
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        result *= a(i, i);
    }
    return result;
}

template <typename T>
MatrixDiagonal<T, MatrixMxN<T>> MatrixMxN<T>::diagonal() const {
    return MatrixDiagonal<T, MatrixMxN>(*this, true);
}

template <typename T>
MatrixDiagonal<T, MatrixMxN<T>> MatrixMxN<T>::offDiagonal() const {
    return MatrixDiagonal<T, MatrixMxN>(*this, false);
}

template <typename T>
MatrixTriangular<T, MatrixMxN<T>> MatrixMxN<T>::strictLowerTri() const {
    return MatrixTriangular<T, MatrixMxN<T>>(*this, false, true);
}

template <typename T>
MatrixTriangular<T, MatrixMxN<T>> MatrixMxN<T>::strictUpperTri() const {
    return MatrixTriangular<T, MatrixMxN<T>>(*this, true, true);
}

template <typename T>
MatrixTriangular<T, MatrixMxN<T>> MatrixMxN<T>::lowerTri() const {
    return MatrixTriangular<T, MatrixMxN<T>>(*this, false, false);
}

template <typename T>
MatrixTriangular<T, MatrixMxN<T>> MatrixMxN<T>::upperTri() const {
    return MatrixTriangular<T, MatrixMxN<T>>(*this, true, false);
}

template <typename T>
MatrixMxN<T> MatrixMxN<T>::transposed() const {
    MatrixMxN mt(cols(), rows());
    parallelForEachIndex([&](size_t i, size_t j) { mt(j, i) = (*this)(i, j); });
    return mt;
}

template <typename T>
MatrixMxN<T> MatrixMxN<T>::inverse() const {
    MatrixMxN mInv(*this);
    mInv.invert();
    return mInv;
}

template <typename T>
template <typename U>
MatrixTypeCast<U, MatrixMxN<T>, T> MatrixMxN<T>::castTo() const {
    return MatrixTypeCast<U, MatrixMxN, T>(*this);
}

template <typename T>
template <typename E>
MatrixMxN<T>& MatrixMxN<T>::operator=(const E& m) {
    set(m);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator=(const MatrixMxN& other) {
    set(other);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator=(MatrixMxN&& other) {
    _elements = std::move(other._elements);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator+=(const T& s) {
    iadd(s);
    return *this;
}

template <typename T>
template <typename E>
MatrixMxN<T>& MatrixMxN<T>::operator+=(const E& m) {
    iadd(m);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator-=(const T& s) {
    isub(s);
    return *this;
}

template <typename T>
template <typename E>
MatrixMxN<T>& MatrixMxN<T>::operator-=(const E& m) {
    isub(m);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator*=(const T& s) {
    imul(s);
    return *this;
}

template <typename T>
template <typename E>
MatrixMxN<T>& MatrixMxN<T>::operator*=(const E& m) {
    imul(m);
    return *this;
}

template <typename T>
MatrixMxN<T>& MatrixMxN<T>::operator/=(const T& s) {
    idiv(s);
    return *this;
}

template <typename T>
T& MatrixMxN<T>::operator[](size_t i) {
    return _elements[i];
}

template <typename T>
const T& MatrixMxN<T>::operator[](size_t i) const {
    return _elements[i];
}

template <typename T>
T& MatrixMxN<T>::operator()(size_t i, size_t j) {
    return _elements(j, i);
}

template <typename T>
const T& MatrixMxN<T>::operator()(size_t i, size_t j) const {
    return _elements(j, i);
}

template <typename T>
template <typename E>
bool MatrixMxN<T>::operator==(const MatrixExpression<T, E>& m) const {
    return isEqual(m);
}

template <typename T>
template <typename E>
bool MatrixMxN<T>::operator!=(const MatrixExpression<T, E>& m) const {
    return !isEqual(m);
}

template <typename T>
template <typename Callback>
void MatrixMxN<T>::forEach(Callback func) const {
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            func((*this)(i, j));
        }
    }
}

template <typename T>
template <typename Callback>
void MatrixMxN<T>::forEachIndex(Callback func) const {
    for (size_t i = 0; i < rows(); ++i) {
        for (size_t j = 0; j < cols(); ++j) {
            func(i, j);
        }
    }
}

template <typename T>
template <typename Callback>
void MatrixMxN<T>::parallelForEach(Callback func) {
    parallelFor(kZeroSize, cols(), kZeroSize, rows(),
                [&](size_t j, size_t i) { func((*this)(i, j)); });
}

template <typename T>
template <typename Callback>
void MatrixMxN<T>::parallelForEachIndex(Callback func) const {
    parallelFor(kZeroSize, cols(), kZeroSize, rows(),
                [&](size_t j, size_t i) { func(i, j); });
}

// MARK: Builders

template <typename T>
MatrixConstant<T> MatrixMxN<T>::makeZero(size_t m, size_t n) {
    return MatrixConstant<T>(m, n, 0);
}

template <typename T>
MatrixIdentity<T> MatrixMxN<T>::makeIdentity(size_t m) {
    return MatrixIdentity<T>(m);
}

}  // namespace jet

#endif  // INCLUDE_JET_DETAIL_MATRIX_MXN_INL_H_
