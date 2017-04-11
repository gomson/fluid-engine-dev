// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_DETAIL_MATRIX_CSR_INL_H_
#define INCLUDE_JET_DETAIL_MATRIX_CSR_INL_H_

#include <jet/cpp_utils.h>
#include <jet/math_utils.h>
#include <jet/matrix_csr.h>
#include <jet/parallel.h>

#include <numeric>

namespace jet {

template <typename T, typename VE>
MatrixCsrVectorMul<T, VE>::MatrixCsrVectorMul(const MatrixCsr<T>& m,
                                              const VE& v)
    : _m(m), _v(v) {
    JET_ASSERT(_m.cols() == _v.size());
}

template <typename T, typename VE>
size_t MatrixCsrVectorMul<T, VE>::size() const {
    return _v.size();
}

template <typename T, typename VE>
T MatrixCsrVectorMul<T, VE>::operator[](size_t i) const {
    auto rp = _m.rowPointersBegin();
    auto ci = _m.columnIndicesBegin();

    size_t colBegin = rp[i];
    size_t colEnd = rp[i + 1];

    T sum = 0;

    for (size_t jj = colBegin; jj < colEnd; ++jj) {
        size_t j = ci[jj];
        sum += rp[j] * _v[j];
    }

    return sum;
}

//

template <typename T>
MatrixCsr<T>::Element::Element() : i(0), j(0), value(0) {}

template <typename T>
MatrixCsr<T>::Element::Element(size_t i_, size_t j_, const T& value_)
    : i(i_), j(j_), value(value_) {}

//

template <typename T>
MatrixCsr<T>::MatrixCsr() {
    _rowPointers.push_back(0);
}

template <typename T>
MatrixCsr<T>::MatrixCsr(
    const std::initializer_list<std::initializer_list<T>>& lst, T epsilon) {
    compress(lst, epsilon);
}

template <typename T>
template <typename E>
MatrixCsr<T>::MatrixCsr(const MatrixExpression<T, E>& other, T epsilon) {
    compress(other, epsilon);
}

template <typename T>
MatrixCsr<T>::MatrixCsr(const MatrixCsr& other) {
    set(other);
}

template <typename T>
MatrixCsr<T>::MatrixCsr(MatrixCsr&& other) {
    (*this) = std::move(other);
}

template <typename T>
void MatrixCsr<T>::set(const T& s) {
    std::fill(_nonZeros.begin(), _nonZeros.end(), s);
}

template <typename T>
void MatrixCsr<T>::set(const MatrixCsr& other) {
    _size = other._size;
    _nonZeros = other._nonZeros;
    _rowPointers = other._rowPointers;
    _columnIndices = other._columnIndices;
}

template <typename T>
void MatrixCsr<T>::compress(
    const std::initializer_list<std::initializer_list<T>>& lst, T epsilon) {
    size_t numRows = lst.size();
    size_t numCols = (numRows > 0) ? lst.begin()->size() : 0;

    _size = {numRows, numCols};
    _nonZeros.clear();
    _rowPointers.resize(numRows + 1);
    _columnIndices.clear();

    auto rowIter = lst.begin();
    for (size_t j = 0; j < numRows; ++j) {
        JET_ASSERT(numCols == rowIter->size());
        _rowPointers.push_back(_nonZeros.size());

        auto colIter = rowIter->begin();
        for (size_t i = 0; i < numCols; ++i) {
            if (std::fabs(*colIter) > epsilon) {
                _nonZeros.push_back(*colIter);
                _columnIndices.push_back(j);
            }

            ++colIter;
        }
        ++rowIter;
    }

    _rowPointers.push_back(_nonZeros.size());
}

template <typename T>
template <typename E>
void MatrixCsr<T>::compress(const MatrixExpression<T, E>& other, T epsilon) {
    size_t numRows = other.rows();
    size_t numCols = other.cols();

    _size = {numRows, numCols};
    _nonZeros.clear();
    _columnIndices.clear();

    const E& expression = other();

    for (size_t i = 0; i < numRows; ++i) {
        _rowPointers.push_back(_nonZeros.size());

        for (size_t j = 0; j < numCols; ++j) {
            T val = expression(i, j);
            if (std::fabs(val) > epsilon) {
                _nonZeros.push_back(val);
                _columnIndices.push_back(j);
            }
        }
    }

    _rowPointers.push_back(_nonZeros.size());
}

template <typename T>
void MatrixCsr<T>::addElement(size_t i, size_t j, const T& value) {
    addElement({i, j, value});
}

template <typename T>
void MatrixCsr<T>::addElement(const Element& element) {
    ssize_t numRowsToAdd = (ssize_t)element.i - (ssize_t)_size.x + 1;
    if (numRowsToAdd > 0) {
        for (size_t i = 0; i < numRowsToAdd; ++i) {
            addRow({}, {});
        }
    }

    _size.y = std::max(_size.y, element.j + 1);

    size_t rowBegin = _rowPointers[element.i];
    size_t rowEnd = _rowPointers[element.i + 1];

    auto colIdxIter =
        std::lower_bound(_columnIndices.begin() + rowBegin,
                         _columnIndices.begin() + rowEnd, element.j);
    auto offset = colIdxIter - _columnIndices.begin();

    _columnIndices.insert(colIdxIter, element.j);
    _nonZeros.insert(_nonZeros.begin() + offset, element.value);

    for (size_t i = element.i + 1; i < _rowPointers.size(); ++i) {
        ++_rowPointers[i];
    }
}

template <typename T>
void MatrixCsr<T>::addRow(const NonZeroContainerType& nonZeros,
                          const IndexContainerType& columnIndices) {
    JET_ASSERT(nonZeros.size() == columnIndices.size());

    ++_size.x;

    // TODO: Implement zip iterator
    std::vector<std::pair<T, size_t>> zipped;
    for (size_t i = 0; i < nonZeros.size(); ++i) {
        zipped.emplace_back(nonZeros[i], columnIndices[i]);
        _size.y = std::max(_size.y, columnIndices[i] + 1);
    }
    std::sort(zipped.begin(), zipped.end(),
              [](std::pair<T, size_t> a, std::pair<T, size_t> b) {
                  return a.second < b.second;
              });
    for (size_t i = 0; i < zipped.size(); ++i) {
        _nonZeros.push_back(zipped[i].first);
        _columnIndices.push_back(zipped[i].second);
    }

    _rowPointers.push_back(_nonZeros.size());
}

template <typename T>
void MatrixCsr<T>::setElement(size_t i, size_t j, const T& value) {
    setElement({i, j, value});
}

template <typename T>
void MatrixCsr<T>::setElement(const Element& element) {
    size_t nzIndex = hasElement(element.i, element.j);
    if (nzIndex == kMaxSize) {
        addElement(element);
    } else {
        _nonZeros[nzIndex] = element.value;
    }
}

template <typename T>
bool MatrixCsr<T>::isEqual(const MatrixCsr& other) const {
    if (_size != other._size) {
        return false;
    }

    if (_nonZeros.size() != other._nonZeros.size()) {
        return false;
    }

    for (size_t i = 0; i < _nonZeros.size(); ++i) {
        if (_nonZeros[i] != other._nonZeros[i]) {
            return false;
        }
        if (_columnIndices[i] != other._columnIndices[i]) {
            return false;
        }
    }

    for (size_t i = 0; i < _rowPointers.size(); ++i) {
        if (_rowPointers[i] != other._rowPointers[i]) {
            return false;
        }
    }

    return true;
}

template <typename T>
bool MatrixCsr<T>::isSimilar(const MatrixCsr& other, double tol) const {
    if (_size != other._size) {
        return false;
    }

    if (_nonZeros.size() != other._nonZeros.size()) {
        return false;
    }

    for (size_t i = 0; i < _nonZeros.size(); ++i) {
        if (std::fabs(_nonZeros[i] - other._nonZeros[i]) > tol) {
            return false;
        }
        if (_columnIndices[i] != other._columnIndices[i]) {
            return false;
        }
    }

    for (size_t i = 0; i < _rowPointers.size(); ++i) {
        if (_rowPointers[i] != other._rowPointers[i]) {
            return false;
        }
    }

    return true;
}

template <typename T>
bool MatrixCsr<T>::isSquare() const {
    return rows() == cols();
}

template <typename T>
Size2 MatrixCsr<T>::size() const {
    return _size;
}

template <typename T>
size_t MatrixCsr<T>::rows() const {
    return _size.x;
}

template <typename T>
size_t MatrixCsr<T>::cols() const {
    return _size.y;
}

template <typename T>
size_t MatrixCsr<T>::numberOfNonZeros() const {
    return _nonZeros.size();
}

template <typename T>
T* MatrixCsr<T>::nonZeroData() {
    return _nonZeros.data();
}

template <typename T>
const T* const MatrixCsr<T>::nonZeroData() const {
    return _nonZeros.data();
}

template <typename T>
size_t* MatrixCsr<T>::rowPointersData() {
    return _rowPointers.data();
}

template <typename T>
const size_t* const MatrixCsr<T>::rowPointersData() const {
    return _rowPointers.data();
}

template <typename T>
size_t* MatrixCsr<T>::columnIndicesData() {
    return _columnIndices.data();
}

template <typename T>
const size_t* const MatrixCsr<T>::columnIndicesData() const {
    return _columnIndices.data();
}

template <typename T>
typename MatrixCsr<T>::ConstNonZeroIterator MatrixCsr<T>::nonZeroBegin() {
    return _nonZeros.begin();
}

template <typename T>
typename MatrixCsr<T>::ConstNonZeroIterator MatrixCsr<T>::nonZeroBegin() const {
    return _nonZeros.cbegin();
}

template <typename T>
typename MatrixCsr<T>::NonZeroIterator MatrixCsr<T>::nonZeroEnd() {
    return _nonZeros.end();
}

template <typename T>
typename MatrixCsr<T>::ConstNonZeroIterator MatrixCsr<T>::nonZeroEnd() const {
    return _nonZeros.cend();
}

template <typename T>
typename MatrixCsr<T>::IndexIterator MatrixCsr<T>::rowPointersBegin() {
    return _rowPointers.begin();
}

template <typename T>
typename MatrixCsr<T>::ConstIndexIterator MatrixCsr<T>::rowPointersBegin()
    const {
    return _rowPointers.cbegin();
}

template <typename T>
typename MatrixCsr<T>::IndexIterator MatrixCsr<T>::rowPointersEnd() {
    return _rowPointers.end();
}

template <typename T>
typename MatrixCsr<T>::ConstIndexIterator MatrixCsr<T>::rowPointersEnd() const {
    return _rowPointers.cend();
}

template <typename T>
typename MatrixCsr<T>::IndexIterator MatrixCsr<T>::columnIndicesBegin() {
    return _columnIndices.begin();
}

template <typename T>
typename MatrixCsr<T>::ConstIndexIterator MatrixCsr<T>::columnIndicesBegin()
    const {
    return _columnIndices.cbegin();
}

template <typename T>
typename MatrixCsr<T>::IndexIterator MatrixCsr<T>::columnIndicesEnd() {
    return _columnIndices.end();
}

template <typename T>
typename MatrixCsr<T>::ConstIndexIterator MatrixCsr<T>::columnIndicesEnd()
    const {
    return _columnIndices.cend();
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::add(const T& s) const {
    MatrixCsr ret(*this);
    parallelFor(kZeroSize, ret._nonZeros.size(),
                [&](size_t i) { ret._nonZeros[i] += s; });
    return ret;
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::add(const MatrixCsr& m) const {
    JET_ASSERT(_size == m._size);

    MatrixCsr ret;
    for (size_t i = 0; i < _size.x; ++i) {
        // TODO: Implement
    }

    return ret;
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::sub(const T& s) const {
    MatrixCsr ret(*this);
    parallelFor(kZeroSize, ret._nonZeros.size(),
                [&](size_t i) { ret._nonZeros[i] -= s; });
    return ret;
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::sub(const MatrixCsr& m) const {
    JET_ASSERT(_size == m._size);

    MatrixCsr ret;
    for (size_t i = 0; i < _size.x; ++i) {
        // TODO: Implement
    }

    return ret;
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::mul(const T& s) const {
    MatrixCsr ret(*this);
    parallelFor(kZeroSize, ret._nonZeros.size(),
                [&](size_t i) { ret._nonZeros[i] *= s; });
    return ret;
}

template <typename T>
template <typename VE>
MatrixCsrVectorMul<T, VE> MatrixCsr<T>::mul(
    const VectorExpression<T, VE>& v) const {
    return MatrixCsrVectorMul<T, VE>(*this, v);
};

template <typename T>
template <typename E>
MatrixCsr<T> MatrixCsr<T>::mul(const E& m) const {
    // TODO: Implement
    (void)m;
    return MatrixCsr();
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::mul(const MatrixCsr& m) const {
    // TODO: Implement
    (void)m;
    return MatrixCsr();
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::div(const T& s) const {
    MatrixCsr ret(*this);
    parallelFor(kZeroSize, ret._nonZeros.size(),
                [&](size_t i) { ret._nonZeros[i] /= s; });
    return ret;
}

template <typename T>
T MatrixCsr<T>::sum() const {
    return parallelReduce(kZeroSize, numberOfNonZeros(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result += _nonZeros[i];
                              }
                              return result;
                          },
                          std::plus<T>());
}

template <typename T>
T MatrixCsr<T>::avg() const {
    return sum() / numberOfNonZeros();
}

template <typename T>
T MatrixCsr<T>::min() const {
    const T& (*_min)(const T&, const T&) = std::min<T>;
    return parallelReduce(kZeroSize, numberOfNonZeros(),
                          std::numeric_limits<T>::max(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = std::min(result, _nonZeros[i]);
                              }
                              return result;
                          },
                          _min);
}

template <typename T>
T MatrixCsr<T>::max() const {
    const T& (*_max)(const T&, const T&) = std::max<T>;
    return parallelReduce(kZeroSize, numberOfNonZeros(),
                          std::numeric_limits<T>::min(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = std::max(result, _nonZeros[i]);
                              }
                              return result;
                          },
                          _max);
}

template <typename T>
T MatrixCsr<T>::absmin() const {
    return parallelReduce(kZeroSize, numberOfNonZeros(),
                          std::numeric_limits<T>::max(),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = jet::absmin(result, _nonZeros[i]);
                              }
                              return result;
                          },
                          jet::absmin<T>);
}

template <typename T>
T MatrixCsr<T>::absmax() const {
    return parallelReduce(kZeroSize, numberOfNonZeros(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result = jet::absmax(result, _nonZeros[i]);
                              }
                              return result;
                          },
                          jet::absmax<T>);
}

template <typename T>
T MatrixCsr<T>::trace() const {
    JET_ASSERT(isSquare());
    return parallelReduce(kZeroSize, rows(), T(0),
                          [&](size_t start, size_t end, T init) {
                              T result = init;
                              for (size_t i = start; i < end; ++i) {
                                  result += (*this)(i, i);
                              }
                              return result;
                          },
                          std::plus<T>());
}

template <typename T>
T MatrixCsr<T>::determinant() const {
    // TODO: Implement;
    return 0;
}

template <typename T>
template <typename E>
MatrixCsr<T>& MatrixCsr<T>::operator=(const E& m) {
    set(m);
}

template <typename T>
MatrixCsr<T>& MatrixCsr<T>::operator=(const MatrixCsr& other) {
    set(other);
}

template <typename T>
MatrixCsr<T>& MatrixCsr<T>::operator=(MatrixCsr&& other) {
    _size = other._size;
    other._size = Size2();
    _nonZeros = std::move(other._nonZeros);
    _rowPointers = std::move(other._rowPointers);
    _columnIndices = std::move(other._columnIndices);
    return *this;
}

template <typename T>
T MatrixCsr<T>::operator()(size_t i, size_t j) const {
    size_t nzIndex = hasElement(i, j);
    if (nzIndex == kMaxSize) {
        return 0.0;
    } else {
        return _nonZeros[nzIndex];
    }
}

template <typename T>
bool MatrixCsr<T>::operator==(const MatrixCsr& m) const {
    return isEqual(m);
}

template <typename T>
bool MatrixCsr<T>::operator!=(const MatrixCsr& m) const {
    return !isEqual(m);
}

template <typename T>
MatrixCsr<T> MatrixCsr<T>::makeIdentity(size_t m) {
    MatrixCsr ret;
    ret._size = Size2(m, m);
    ret._nonZeros.resize(m, 1.0);
    ret._columnIndices.resize(m);
    std::iota(ret._columnIndices.begin(), ret._columnIndices.end(), kZeroSize);
    ret._rowPointers.resize(m + 1);
    std::iota(ret._rowPointers.begin(), ret._rowPointers.end(), kZeroSize);
    return ret;
}

template <typename T>
size_t MatrixCsr<T>::hasElement(size_t i, size_t j) const {
    if (i >= _size.x || j >= _size.y) {
        return kMaxSize;
    }

    size_t rowBegin = _rowPointers[i];
    size_t rowEnd = _rowPointers[i + 1];

    auto iter = binaryFind(_columnIndices.begin() + rowBegin,
                           _columnIndices.begin() + rowEnd, j);
    if (iter != _columnIndices.begin() + rowEnd) {
        return static_cast<size_t>(iter - _columnIndices.begin());
    } else {
        return kMaxSize;
    }
}

}  // namespace jet

#endif  // INCLUDE_JET_DETAIL_MATRIX_CSR_INL_H_
