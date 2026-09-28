// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include <vector>
#include <complex>
#include <mkl.h>
#include "ComplexScalar.hpp"
#include "ComplexVector.hpp"
#include "floatmatrix.h"
#include <fstream>
#include <iomanip>  // 如果需要设置输出精度

class ComplexMatrix
{
private:
    int rows_, cols_;
    std::vector<MKL_Complex16> data_; // 列主序存储

public:
    // 构造函数
    ComplexMatrix(int rows = 0, int cols = 0);
    ComplexMatrix(int rows, int cols, const ComplexScalar &init_val);
    ComplexMatrix(int rows, int cols, FloatArray &real, FloatArray &imag);
    ComplexMatrix(const FloatMatrix &realPart, const FloatMatrix &imagPart);
    ComplexMatrix(const FloatMatrix &data,std::string type);

    // 获取尺寸
    int Rows() const { return rows_; }
    int Cols() const { return cols_; }
    int Size() const { return rows_ * cols_; }

    // 访问元素（列主序）
    MKL_Complex16 &operator()(int i, int j);
    const MKL_Complex16 &operator()(int i, int j) const;

    // 返回 ComplexScalar（值类型）
    ComplexScalar At(int i, int j) const;

    // 调整矩阵大小并初始化（可指定初值，默认零）
    void resize(int rows, int cols, const ComplexScalar &init_val = ComplexScalar(0.0, 0.0));

    ComplexMatrix ConcatenateRow(const ComplexMatrix& A, const ComplexMatrix& B);
    ComplexMatrix ConcatenateCol(const ComplexMatrix& A, const ComplexMatrix& B);
    ComplexMatrix ConcatenateDiag(const ComplexMatrix& A, const ComplexMatrix& B);

    void EliminateElement(const IntArray& index);

    ComplexMatrix transpose();

    ComplexVector GetColumn(int col) const;

    ComplexMatrix GetRow(const IntArray& row) const;

    void SetSubMatrix(const IntArray& row,const IntArray& col, const ComplexMatrix& sub);
    ComplexMatrix SetColumn(int col, const ComplexVector &value);

    IntArray FindMaxElementLoaction2();


    // 初始化为零
    void Zero();

    // 返回非零元素的位置
    void findNonZero(IntArray& rows, IntArray& cols) const;

    // 复制构造与赋值
    ComplexMatrix(const ComplexMatrix &other);
    ComplexMatrix &operator=(const ComplexMatrix &other);

    // 设置指定位置的值
    void Set(int i, int j, const ComplexScalar &value)
    {
        if(i < 0 || i >= rows_ || j < 0 || j >= cols_)
            throw std::out_of_range("ComplexMatrix::Set index out of range");
        data_[i + j * rows_] = { value.real(), value.imag() };
    }

    // 也可以提供直接设置实部和虚部的方法
    void Set(int i, int j, double real, double imag)
    {
        if(i < 0 || i >= rows_ || j < 0 || j >= cols_)
            throw std::out_of_range("ComplexMatrix::Set index out of range");
        data_[i + j * rows_] = { real, imag };
    }



    // 运算符重载
    ComplexMatrix operator+(const ComplexMatrix &other) const;
    ComplexMatrix operator-(const ComplexMatrix &other) const;
    ComplexMatrix operator*(const ComplexMatrix &other) const;
    ComplexMatrix operator*(const ComplexScalar &alpha) const;

    // 矩阵向量乘
    ComplexVector operator*(const ComplexVector &vec) const;

    void ComputeEigenValueAndEigenVector(ComplexVector& EigenValue,ComplexMatrix& EigenVector);

    // 输出
    void Print(std::ostream &out = std::cout) const;

    // MKL 内部数据访问
    MKL_Complex16 *Data() { return data_.data(); }
    const MKL_Complex16 *Data() const { return data_.data(); }

    friend std::ostream &operator<<(std::ostream &out, const ComplexMatrix &mat)
    {
        std::vector<int> col_width(mat.cols_, 0);

        // lambda 生成高精度复数字符串
        auto format_complex = [](const MKL_Complex16 &v) {
            std::ostringstream oss;

            if (v.imag >= 0)
                oss << v.real << "+" << v.imag << "i";
            else
                oss << v.real << v.imag << "i";
            return oss.str();
        };

        // 计算每列最大宽度
        for (int j = 0; j < mat.cols_; ++j)
        {
            for (int i = 0; i < mat.rows_; ++i)
            {
                int len = static_cast<int>(format_complex(mat(i, j)).length());
                if (len > col_width[j])
                    col_width[j] = len;
            }
        }

        // 输出矩阵
        for (int i = 0; i < mat.rows_; ++i)
        {
            for (int j = 0; j < mat.cols_; ++j)
            {
                out << std::setw(col_width[j]) << format_complex(mat(i, j));
                if (j != mat.cols_ - 1)
                    out << "  "; // 列间隔
            }
            out << std::endl;
        }

        return out;
    }


};
