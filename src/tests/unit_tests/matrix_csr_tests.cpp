// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <jet/matrix_csr.h>
#include <jet/matrix_mxn.h>

#include <gtest/gtest.h>

using namespace jet;

TEST(MatrixCsr, Constructors) {
    const MatrixCsrD emptyMat;

    EXPECT_EQ(0u, emptyMat.rows());
    EXPECT_EQ(0u, emptyMat.cols());
    EXPECT_EQ(0u, emptyMat.numberOfNonZeros());
    EXPECT_EQ(emptyMat.nonZeroBegin(), emptyMat.nonZeroEnd());
    EXPECT_EQ(1, emptyMat.rowPointersEnd() - emptyMat.rowPointersBegin());
    EXPECT_EQ(emptyMat.columnIndicesBegin(), emptyMat.columnIndicesEnd());

    const MatrixCsrD matInitLst = {
        {1.0, 0.0, 0.0, -3.0}, {0.0, 3.0, -5.0, 1.0}, {-4.0, 0.0, 1.0, 5.0}};
    EXPECT_EQ(3u, matInitLst.rows());
    EXPECT_EQ(4u, matInitLst.cols());
    EXPECT_EQ(8u, matInitLst.numberOfNonZeros());

    auto iterInitLst = matInitLst.nonZeroBegin();
    EXPECT_EQ(1.0, iterInitLst[0]);
    EXPECT_EQ(-3.0, iterInitLst[1]);
    EXPECT_EQ(3.0, iterInitLst[2]);
    EXPECT_EQ(-5.0, iterInitLst[3]);
    EXPECT_EQ(1.0, iterInitLst[4]);
    EXPECT_EQ(-4.0, iterInitLst[5]);
    EXPECT_EQ(1.0, iterInitLst[6]);
    EXPECT_EQ(5.0, iterInitLst[7]);

    const MatrixMxND matDense = {
        {1.0, 0.01, 0.0, -3.0}, {0.01, 3.0, -5.0, 1.0}, {-4.0, 0.01, 1.0, 5.0}};
    const MatrixCsrD matSparse(matDense, 0.02);
    EXPECT_EQ(3u, matSparse.rows());
    EXPECT_EQ(4u, matSparse.cols());
    EXPECT_EQ(8u, matSparse.numberOfNonZeros());

    auto iterSparse = matSparse.nonZeroBegin();
    EXPECT_EQ(1.0, iterSparse[0]);
    EXPECT_EQ(-3.0, iterSparse[1]);
    EXPECT_EQ(3.0, iterSparse[2]);
    EXPECT_EQ(-5.0, iterSparse[3]);
    EXPECT_EQ(1.0, iterSparse[4]);
    EXPECT_EQ(-4.0, iterSparse[5]);
    EXPECT_EQ(1.0, iterSparse[6]);
    EXPECT_EQ(5.0, iterSparse[7]);

    MatrixCsrD matCopied = matSparse;
    EXPECT_EQ(3u, matCopied.rows());
    EXPECT_EQ(4u, matCopied.cols());
    EXPECT_EQ(8u, matCopied.numberOfNonZeros());

    auto iterCopied = matCopied.nonZeroBegin();
    EXPECT_EQ(1.0, iterCopied[0]);
    EXPECT_EQ(-3.0, iterCopied[1]);
    EXPECT_EQ(3.0, iterCopied[2]);
    EXPECT_EQ(-5.0, iterCopied[3]);
    EXPECT_EQ(1.0, iterCopied[4]);
    EXPECT_EQ(-4.0, iterCopied[5]);
    EXPECT_EQ(1.0, iterCopied[6]);
    EXPECT_EQ(5.0, iterCopied[7]);

    const MatrixCsrD matMoved = std::move(matCopied);
    EXPECT_EQ(3u, matMoved.rows());
    EXPECT_EQ(4u, matMoved.cols());
    EXPECT_EQ(8u, matMoved.numberOfNonZeros());

    auto iterMovied = matMoved.nonZeroBegin();
    EXPECT_EQ(1.0, iterMovied[0]);
    EXPECT_EQ(-3.0, iterMovied[1]);
    EXPECT_EQ(3.0, iterMovied[2]);
    EXPECT_EQ(-5.0, iterMovied[3]);
    EXPECT_EQ(1.0, iterMovied[4]);
    EXPECT_EQ(-4.0, iterMovied[5]);
    EXPECT_EQ(1.0, iterMovied[6]);
    EXPECT_EQ(5.0, iterMovied[7]);

    EXPECT_EQ(0u, matCopied.rows());
    EXPECT_EQ(0u, matCopied.cols());
    EXPECT_EQ(0u, matCopied.numberOfNonZeros());
    EXPECT_EQ(matCopied.nonZeroBegin(), matCopied.nonZeroEnd());
    EXPECT_EQ(matCopied.rowPointersBegin(), matCopied.rowPointersEnd());
    EXPECT_EQ(matCopied.columnIndicesBegin(), matCopied.columnIndicesEnd());
}

TEST(MatrixCsr, BasicSetters) {
    // Compress initializer list
    const std::initializer_list<std::initializer_list<double>> initLst = {
        {1.0, 0.01, 0.0, -3.0}, {0.01, 3.0, -5.0, 1.0}, {-4.0, 0.01, 1.0, 5.0}};
    MatrixCsrD matInitLst;
    matInitLst.compress(initLst, 0.02);
    EXPECT_EQ(3u, matInitLst.rows());
    EXPECT_EQ(4u, matInitLst.cols());
    EXPECT_EQ(8u, matInitLst.numberOfNonZeros());

    auto iterInitLst = matInitLst.nonZeroBegin();
    EXPECT_EQ(1.0, iterInitLst[0]);
    EXPECT_EQ(-3.0, iterInitLst[1]);
    EXPECT_EQ(3.0, iterInitLst[2]);
    EXPECT_EQ(-5.0, iterInitLst[3]);
    EXPECT_EQ(1.0, iterInitLst[4]);
    EXPECT_EQ(-4.0, iterInitLst[5]);
    EXPECT_EQ(1.0, iterInitLst[6]);
    EXPECT_EQ(5.0, iterInitLst[7]);

    // Set scalar
    matInitLst.set(42.0);
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(42.0, iterInitLst[i]);
    }

    // Compress dense matrix
    const MatrixMxND matDense = {
        {1.0, 0.01, 0.0, -3.0}, {0.01, 3.0, -5.0, 1.0}, {-4.0, 0.01, 1.0, 5.0}};
    MatrixCsrD matSparse;
    matSparse.compress(matDense, 0.02);
    EXPECT_EQ(3u, matSparse.rows());
    EXPECT_EQ(4u, matSparse.cols());
    EXPECT_EQ(8u, matSparse.numberOfNonZeros());

    auto iterSparse = matSparse.nonZeroBegin();
    EXPECT_EQ(1.0, iterSparse[0]);
    EXPECT_EQ(-3.0, iterSparse[1]);
    EXPECT_EQ(3.0, iterSparse[2]);
    EXPECT_EQ(-5.0, iterSparse[3]);
    EXPECT_EQ(1.0, iterSparse[4]);
    EXPECT_EQ(-4.0, iterSparse[5]);
    EXPECT_EQ(1.0, iterSparse[6]);
    EXPECT_EQ(5.0, iterSparse[7]);

    // Set other CSR mat
    matInitLst.set(matSparse);
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(iterSparse[i], iterInitLst[i]);
    }

    // Add/set element
    MatrixCsrD matAddElem;
    matAddElem.addElement(0, 0, 1.0);
    matAddElem.setElement(0, 3, -3.0);
    matAddElem.addElement(1, 1, 3.0);
    matAddElem.setElement(1, 2, -5.0);
    matAddElem.addElement(1, 3, 1.0);
    matAddElem.setElement(2, 0, -4.0);
    matAddElem.addElement(2, 2, 1.0);
    matAddElem.setElement(2, 3, 5.0);

    EXPECT_EQ(3u, matAddElem.rows());
    EXPECT_EQ(4u, matAddElem.cols());
    EXPECT_EQ(8u, matAddElem.numberOfNonZeros());

    auto iterAddElem = matAddElem.nonZeroBegin();
    EXPECT_EQ(1.0, iterAddElem[0]);
    EXPECT_EQ(-3.0, iterAddElem[1]);
    EXPECT_EQ(3.0, iterAddElem[2]);
    EXPECT_EQ(-5.0, iterAddElem[3]);
    EXPECT_EQ(1.0, iterAddElem[4]);
    EXPECT_EQ(-4.0, iterAddElem[5]);
    EXPECT_EQ(1.0, iterAddElem[6]);
    EXPECT_EQ(5.0, iterAddElem[7]);

    matAddElem.setElement(1, 3, 7.0);
    EXPECT_EQ(7.0, iterAddElem[4]);

    // Add element in random order
    MatrixCsrD matAddElemRandom;
    matAddElemRandom.addElement(2, 2, 1.0);
    matAddElemRandom.addElement(0, 3, -3.0);
    matAddElemRandom.addElement(2, 0, -4.0);
    matAddElemRandom.addElement(1, 1, 3.0);
    matAddElemRandom.addElement(2, 3, 5.0);
    matAddElemRandom.addElement(1, 3, 1.0);
    matAddElemRandom.addElement(1, 2, -5.0);
    matAddElemRandom.addElement(0, 0, 1.0);

    EXPECT_EQ(3u, matAddElemRandom.rows());
    EXPECT_EQ(4u, matAddElemRandom.cols());
    EXPECT_EQ(8u, matAddElemRandom.numberOfNonZeros());

    auto iterAddElemRandom = matAddElemRandom.nonZeroBegin();
    EXPECT_EQ(1.0, iterAddElemRandom[0]);
    EXPECT_EQ(-3.0, iterAddElemRandom[1]);
    EXPECT_EQ(3.0, iterAddElemRandom[2]);
    EXPECT_EQ(-5.0, iterAddElemRandom[3]);
    EXPECT_EQ(1.0, iterAddElemRandom[4]);
    EXPECT_EQ(-4.0, iterAddElemRandom[5]);
    EXPECT_EQ(1.0, iterAddElemRandom[6]);
    EXPECT_EQ(5.0, iterAddElemRandom[7]);

    // Add row
    MatrixCsrD matAddRow;
    matAddRow.addRow({1.0, -3.0}, {0, 3});
    matAddRow.addRow({3.0, -5.0, 1.0}, {1, 2, 3});
    matAddRow.addRow({-4.0, 1.0, 5.0}, {0, 2, 3});

    EXPECT_EQ(3u, matAddRow.rows());
    EXPECT_EQ(4u, matAddRow.cols());
    EXPECT_EQ(8u, matAddRow.numberOfNonZeros());

    auto iterAddRow = matAddRow.nonZeroBegin();
    EXPECT_EQ(1.0, iterAddRow[0]);
    EXPECT_EQ(-3.0, iterAddRow[1]);
    EXPECT_EQ(3.0, iterAddRow[2]);
    EXPECT_EQ(-5.0, iterAddRow[3]);
    EXPECT_EQ(1.0, iterAddRow[4]);
    EXPECT_EQ(-4.0, iterAddRow[5]);
    EXPECT_EQ(1.0, iterAddRow[6]);
    EXPECT_EQ(5.0, iterAddRow[7]);
}

TEST(MatrixCsr, GetterOperators) {
    const MatrixMxND matDense = {
            {1.0, 0.0, 0.0, -3.0}, {0.0, 3.0, -5.0, 1.0}, {-4.0, 0.0, 1.0, 5.0}};
    const MatrixCsrD matSparse(matDense, 0.02);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_EQ(matDense(i, j), matSparse(i, j));
        }
    }
}

TEST(MatrixCsr, Builders) {
    const MatrixCsrD matIden = MatrixCsrD::makeIdentity(5);
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            if (i == j) {
                EXPECT_EQ(1.0, matIden(i, j));
            } else {
                EXPECT_EQ(0.0, matIden(i, j));
            }
        }
    }
}
