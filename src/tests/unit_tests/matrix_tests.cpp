// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <jet/matrix.h>
#include <jet/vector.h>

#include <gtest/gtest.h>

using namespace jet;

TEST(Matrix, Constructors) {
    Matrix<double, 2, 3> mat;

    EXPECT_EQ(2u, mat.rows());
    EXPECT_EQ(3u, mat.cols());

    for (double elem : mat) {
        EXPECT_DOUBLE_EQ(0.0, elem);
    }

    Matrix<double, 2, 3> mat2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);

    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(i + 1.0, mat2[i]);
    }

    Matrix<double, 2, 3> mat3 = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(i + 1.0, mat3[i]);
    }

    Matrix<double, 2, 3> mat4(mat3);

    for (int i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(i + 1.0, mat4[i]);
    }
}

TEST(Matrix, BasicSetters) {
    Matrix<double, 4, 2> mat;
    mat.set(5.0);
    EXPECT_EQ(4u, mat.rows());
    EXPECT_EQ(2u, mat.cols());
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(5.0, mat[i]);
    }

    mat.set(7.0);
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(7.0, mat[i]);
    }

    mat.set({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}});
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(i + 1.0, mat[i]);
    }

    Matrix<double, 4, 2> mat2;
    mat2.set(mat);
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(i + 1.0, mat2[i]);
    }

    mat.setDiagonal(10.0);
    for (size_t i = 0; i < 8; ++i) {
        if (i == 0 || i == 3) {
            EXPECT_EQ(10.0, mat[i]);
        } else {
            EXPECT_EQ(mat2[i], mat[i]);
        }
    }

    mat.setOffDiagonal(-1.0);
    for (size_t i = 0; i < 8; ++i) {
        if (i == 0 || i == 3) {
            EXPECT_EQ(10.0, mat[i]);
        } else {
            EXPECT_EQ(-1.0, mat[i]);
        }
    }

    Vector<double, 2> row = {2.0, 3.0};
    mat.setRow(2, row);
    for (size_t i = 0; i < 8; ++i) {
        if (i == 0 || i == 3) {
            EXPECT_EQ(10.0, mat[i]);
        } else if (i == 4) {
            EXPECT_EQ(2.0, mat[i]);
        } else if (i == 5) {
            EXPECT_EQ(3.0, mat[i]);
        } else {
            EXPECT_EQ(-1.0, mat[i]);
        }
    }

    mat.set({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}});
    mat2.set({{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}});
    EXPECT_TRUE(mat.isEqual(mat2));

    mat2.set({{1.01, 2.01}, {3.01, 4.01}, {4.99, 5.99}, {6.99, 7.99}});
    EXPECT_TRUE(mat.isSimilar(mat2, 0.02));
    EXPECT_FALSE(mat.isSimilar(mat2, 0.005));

    EXPECT_FALSE(mat.isSquare());
}

TEST(Matrix, BinaryOperatorMethod) {
    const Matrix<double, 2, 3> matA = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

    Matrix<double, 2, 3> matB = matA.add(3.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i + 4.5, matB[i]);
    }

    Matrix<double, 2, 3> matC = {{3.0, -1.0, 2.0}, {9.0, 2.0, 8.0}};
    matB = matA.add(matC);
    Matrix<double, 2, 3> ans = {{4.0, 1.0, 5.0}, {13.0, 7.0, 14.0}};
    EXPECT_TRUE(ans.isEqual(matB));

    matB = matA.sub(1.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i - 0.5, matB[i]);
    }

    matB = matA.sub(matC);
    ans = {{-2.0, 3.0, 1.0}, {-5.0, 3.0, -2.0}};
    EXPECT_TRUE(ans.isEqual(matB));

    matB = matA.mul(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(2.0 * (i + 1.0), matB[i]);
    }

    Matrix<double, 3, 2> matD = {{3.0, -1.0}, {2.0, 9.0}, {2.0, 8.0}};
    Matrix<double, 2, 2> matE = matA.mul(matD);
    Matrix<double, 2, 2> ans2 = {{13.0, 41.0}, {34.0, 89.0}};
    EXPECT_TRUE(ans2.isEqual(matE));

    matB = matA.div(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ((i + 1.0) / 2.0, matB[i]);
    }

    matB = matA.radd(3.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i + 4.5, matB[i]);
    }

    matB = matA.rsub(1.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(0.5 - i, matB[i]);
    }

    matC = {{3.0, -1.0, 2.0}, {9.0, 2.0, 8.0}};
    matB = matA.rsub(matC);
    ans = {{2.0, -3.0, -1.0}, {5.0, -3.0, 2.0}};
    EXPECT_EQ(ans, matB);

    matB = matA.rmul(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(2.0 * (i + 1.0), matB[i]);
    }

    matD = {{3.0, -1.0}, {2.0, 9.0}, {2.0, 8.0}};
    matE = matD.rmul(matA);
    ans2 = {{13.0, 41.0}, {34.0, 89.0}};
    EXPECT_EQ(ans2, matE);

    matB = matA.rdiv(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(2.0 / (i + 1.0), matB[i]);
    }
}

TEST(Matrix, AugmentedOperatorMethod) {
    const Matrix<double, 2, 3> matA = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    const Matrix<double, 2, 3> matB = {{3.0, -1.0, 2.0}, {9.0, 2.0, 8.0}};

    Matrix<double, 2, 3> mat = matA;
    mat.iadd(3.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i + 4.5, mat[i]);
    }

    mat = matA;
    mat.iadd(matB);
    Matrix<double, 2, 3> ans = {{4.0, 1.0, 5.0}, {13.0, 7.0, 14.0}};
    EXPECT_EQ(ans, mat);

    mat = matA;
    mat.isub(1.5);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i - 0.5, mat[i]) << i;
    }

    mat = matA;
    mat.isub(matB);
    ans = {{-2.0, 3.0, 1.0}, {-5.0, 3.0, -2.0}};
    EXPECT_EQ(ans, mat);

    mat = matA;
    mat.imul(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(2.0 * (i + 1.0), mat[i]);
    }

    Matrix<double, 2, 2> matA2 = {{1.0, 2.0}, {4.0, 5.0}};
    const Matrix<double, 2, 2> matC2 = {{3.0, -1.0}, {2.0, 9.0}};
    matA2.imul(matC2);

    const Matrix<double, 2, 2> ans2 = {{7.0, 17.0}, {22.0, 41.0}};
    EXPECT_EQ(ans2, matA2);

    mat = matA;
    mat.idiv(2.0);
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ((i + 1.0) / 2.0, mat[i]);
    }
}

TEST(Matrix, ComplexGetters) {
    const Matrix<double, 2, 3> matA = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

    EXPECT_EQ(21.0, matA.sum());
    EXPECT_DOUBLE_EQ(21.0 / 6.0, matA.avg());

    const Matrix<double, 2, 3> matB = {{3.0, -1.0, 2.0}, {-9.0, 2.0, 8.0}};
    EXPECT_EQ(-9.0, matB.min());
    EXPECT_EQ(8.0, matB.max());
    EXPECT_EQ(-1.0, matB.absmin());
    EXPECT_EQ(-9.0, matB.absmax());

    const Matrix<double, 3, 3> matC = {
            {3.0, -1.0, 2.0}, {-9.0, 2.0, 8.0}, {4.0, 3.0, 6.0}};
    EXPECT_EQ(11.0, matC.trace());

    EXPECT_DOUBLE_EQ(-192.0, matC.determinant());

    Matrix<double, 2, 3> mat = matA.diagonal();
    Matrix<double, 2, 3> ans = {{1.0, 0.0, 0.0}, {0.0, 5.0, 0.0}};
    EXPECT_EQ(ans, mat);

    mat = matA.offDiagonal();
    ans = {{0.0, 2.0, 3.0}, {4.0, 0.0, 6.0}};
    EXPECT_EQ(ans, mat);

    mat = matC.strictLowerTri();
    ans = {{0.0, 0.0, 0.0}, {-9.0, 0.0, 0.0}, {4.0, 3.0, 0.0}};
    EXPECT_EQ(ans, mat);

    mat = matC.strictUpperTri();
    ans = {{0.0, -1.0, 2.0}, {0.0, 0.0, 8.0}, {0.0, 0.0, 0.0}};
    EXPECT_EQ(ans, mat);

    mat = matC.lowerTri();
    ans = {{3.0, 0.0, 0.0}, {-9.0, 2.0, 0.0}, {4.0, 3.0, 6.0}};
    EXPECT_EQ(ans, mat);

    mat = matC.upperTri();
    ans = {{3.0, -1.0, 2.0}, {0.0, 2.0, 8.0}, {0.0, 0.0, 6.0}};
    EXPECT_EQ(ans, mat);

    const Matrix<float, 3, 3> matF = matC.castTo<float>();
    const Matrix<float, 3, 3> ansF = {{3.f, -1.f, 2.f}, {-9.f, 2.f, 8.f}, {4.f, 3.f, 6.f}};
    EXPECT_EQ(ansF, matF);

    const Matrix<double, 3, 2> matT = matA.transposed();
    const Matrix<double, 3, 2> ansT = {{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}};
    EXPECT_EQ(ansT, matT);

    Matrix<double, 3, 3> matI = {{1.0, 2.0, 3.0}, {2.0, 5.0, 3.0}, {1.0, 0.0, 8.0}};
    Matrix<double, 3, 3> mat2I = matI.inverse();
    Matrix<double, 3, 3> ansI = {{-40.0, 16.0, 9.0}, {13.0, -5.0, -3.0}, {5.0, -2.0, -1.0}};
    EXPECT_TRUE(mat2I.isSimilar(ansI, 1e-9));

    matI = {{1.0, 2.0, 3.0}, {0.0, 1.0, 4.0}, {5.0, 6.0, 0.0}};
    mat2I = matI.inverse();
    ansI = {{-24.0, 18.0, 5.0}, {20.0, -15.0, -4.0}, {-5.0, 4.0, 1.0}};
    EXPECT_TRUE(mat2I.isSimilar(ansI, 1e-9));

    matI = {{0.0, 1.0, 4.0}, {1.0, 2.0, 3.0}, {5.0, 6.0, 0.0}};
    mat2I = matI.inverse();
    ansI = {{18.0, -24.0, 5.0}, {-15.0, 20.0, -4.0}, {4.0, -5.0, 1.0}};
    EXPECT_TRUE(mat2I.isSimilar(ansI, 1e-9));
}

TEST(Matrix, Modifiers) {
    Matrix<double, 3, 3> mat = {{9.0, -8.0, 7.0}, {-6.0, 5.0, -4.0}, {3.0, -2.0, 1.0}};
    mat.transpose();

    Matrix<double, 3, 3> ans = {{9.0, -6.0, 3.0}, {-8.0, 5.0, -2.0}, {7.0, -4.0, 1.0}};
    EXPECT_EQ(ans, mat);

    mat = {{1.0, 2.0, 3.0}, {2.0, 5.0, 3.0}, {1.0, 0.0, 8.0}};
    mat.invert();
    ans = {{-40.0, 16.0, 9.0}, {13.0, -5.0, -3.0}, {5.0, -2.0, -1.0}};
    EXPECT_TRUE(mat.isSimilar(ans, 1e-9));
}

TEST(Matrix, SetterOperators) {
    const Matrix<double, 2, 3> matA = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    const Matrix<double, 2, 3> matB = {{3.0, -1.0, 2.0}, {9.0, 2.0, 8.0}};

    Matrix<double, 2, 3> mat = matA;
    mat += 3.5;
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i + 4.5, mat[i]);
    }

    mat = matA;
    mat += matB;
    Matrix<double, 2, 3> ans = {{4.0, 1.0, 5.0}, {13.0, 7.0, 14.0}};
    EXPECT_EQ(ans, mat);

    mat = matA;
    mat -= 1.5;
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(i - 0.5, mat[i]) << i;
    }

    mat = matA;
    mat -= matB;
    ans = {{-2.0, 3.0, 1.0}, {-5.0, 3.0, -2.0}};
    EXPECT_EQ(ans, mat);

    mat = matA;
    mat *= 2.0;
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(2.0 * (i + 1.0), mat[i]);
    }

    Matrix<double, 2, 2> mat2 = {{1.0, 2.0}, {4.0, 5.0}};
    const Matrix<double, 2, 2> matC2 = {{3.0, -1.0}, {2.0, 9.0}};
    mat2 *= matC2;
    const Matrix<double, 2, 2> ans2 = {{7.0, 17.0}, {22.0, 41.0}};
    EXPECT_EQ(ans2, mat2);

    mat = matA;
    mat /= 2.0;
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ((i + 1.0) / 2.0, mat[i]);
    }
}

TEST(Matrix, GetterOperator) {
    Matrix<double, 2, 4> mat, mat2;
    mat.set({{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}});
    double cnt = 1.0;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_EQ(cnt, mat(i, j));
            cnt += 1.0;
        }
    }

    for (size_t i = 0; i < 8; ++i) {
        EXPECT_EQ(i + 1.0, mat[i]);
    }

    mat.set({{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}});
    mat2.set({{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}});
    EXPECT_EQ(mat, mat2);
}

TEST(Matrix, Builders) {
    const Matrix<double, 3, 4> mat = Matrix<double, 3, 4>::makeZero();
    for (size_t i = 0; i < 12; ++i) {
        EXPECT_EQ(0.0, mat[i]);
    }

    const Matrix<double, 5, 5> mat2 = Matrix<double, 5, 5>::makeIdentity();
    for (size_t i = 0; i < 25; ++i) {
        if (i % 6 == 0) {
            EXPECT_EQ(1.0, mat2[i]);
        } else {
            EXPECT_EQ(0.0, mat2[i]);
        }
    }
}
