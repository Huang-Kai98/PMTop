// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "ComplexMatrix.h"
#include <iostream>
#include <stdexcept>
#include <cstring>  // for memcpy
#include <algorithm>
#include "floatmatrix.h"

// ======================== 构造函数 ========================
ComplexMatrix::ComplexMatrix(int rows, int cols)
    : rows_(rows), cols_(cols), data_(rows * cols)
{
    Zero();
}

ComplexMatrix::ComplexMatrix(int rows, int cols, const ComplexScalar &init_val)
    : rows_(rows), cols_(cols), data_(rows * cols)
{
    MKL_Complex16 val;
    val.real = init_val.real();
    val.imag = init_val.imag();
    std::fill(data_.begin(), data_.end(), val);
}

ComplexMatrix::ComplexMatrix(int rows, int cols,
                             FloatArray &real,
                             FloatArray &imag)
    : rows_(rows), cols_(cols), data_(rows * cols)
{
    if(real.giveSize() != rows*cols || imag.giveSize() != rows*cols)
        throw std::runtime_error("Real/Imag vectors size mismatch with matrix size");
    for(size_t i = 0; i < data_.size(); i++)
    {
        data_[i].real = real[i];
        data_[i].imag = imag[i];
    }
}

ComplexMatrix::ComplexMatrix(const FloatMatrix &realPart, const FloatMatrix &imagPart) {
    if(realPart.giveNumberOfRows() != imagPart.giveNumberOfRows() ||
           realPart.giveNumberOfColumns() != imagPart.giveNumberOfColumns())
    {
        throw std::invalid_argument("Real and imaginary parts must have same dimensions.");
    }

    rows_ = realPart.giveNumberOfRows();
    cols_ = realPart.giveNumberOfColumns();
    data_.resize(rows_ * cols_);

    // 列主序存储
    for(int j = 0; j < cols_; ++j)
    {
        for(int i = 0; i < rows_; ++i)
        {
            int idx = j * rows_ + i; // 列主序索引
            data_[idx].real = realPart(i, j);
            data_[idx].imag = imagPart(i, j);
        }
    }
}

ComplexMatrix::ComplexMatrix(const FloatMatrix &data,std::string type) {

    rows_ = data.giveNumberOfRows();
    cols_ = data.giveNumberOfColumns();
    size_t n = rows_ * cols_;
    data_.resize(n);

    auto* dest = data_.data();
    const double* src = data.givePointer();


    if (type=="Real") {
        // 编译器可以轻易将此循环自动向量化 (Auto-Vectorization)
        // 甚至可以使用 #pragma omp simd
        for (size_t k = 0; k < n; ++k) {
            dest[k].real = src[k];
            dest[k].imag = 0.0;
        }
    } else {
        for (size_t k = 0; k < n; ++k) {
            dest[k].real = 0.0;
            dest[k].imag = src[k];
        }
    }
    // 列主序存储
    // for(int j = 0; j < cols_; ++j)
    // {
    //     for(int i = 0; i < rows_; ++i)
    //     {
    //         if (type=="Real") {
    //             int idx = j * rows_ + i; // 列主序索引
    //             data_[idx].real = data(i, j);
    //             data_[idx].imag = 0.0;
    //         }
    //         else if (type=="Imag") {
    //             int idx = j * rows_ + i; // 列主序索引
    //             data_[idx].real = 0.0;
    //             data_[idx].imag = data(i, j);
    //         }
    //         else {
    //             std::cerr << "Invalid Matrix type: " << type << "\n";
    //         }
    //     }
    // }
}

// ======================== 元素访问 ========================
MKL_Complex16 &ComplexMatrix::operator()(int i, int j)
{
    if(i < 0 || i >= rows_ || j < 0 || j >= cols_)
        throw std::out_of_range("Matrix index out of range");
    return data_[j * rows_ + i]; // 列主序
}

const MKL_Complex16 &ComplexMatrix::operator()(int i, int j) const
{
    if(i < 0 || i >= rows_ || j < 0 || j >= cols_)
        throw std::out_of_range("Matrix index out of range");
    return data_[j * rows_ + i]; // 列主序
}

// 返回 ComplexScalar（值类型）
ComplexScalar ComplexMatrix::At(int i, int j) const
{
    if(i < 0 || i >= rows_ || j < 0 || j >= cols_)
        throw std::out_of_range("Matrix index out of range");

    const MKL_Complex16 &v = data_[j * rows_ + i]; // 列主序
    return ComplexScalar(v.real, v.imag);
}

void ComplexMatrix::resize(int rows, int cols, const ComplexScalar &init_val)
{
    rows_ = rows;
    cols_ = cols;
    data_.resize(rows * cols);

    MKL_Complex16 val;
    val.real = init_val.real();
    val.imag = init_val.imag();
    std::fill(data_.begin(), data_.end(), val);
}

void ComplexMatrix::EliminateElement(const IntArray &index) {
    for (int i=0; i<index.giveSize(); ++i) {
        for (int j=0; j<cols_; ++j) {
            this->Set(index(i),j,0.0,0.0);
        }
    }
}


ComplexMatrix ComplexMatrix::transpose() {
    ComplexMatrix res(cols_,rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            res(j, i) = (*this)(i, j);
        }
    }
    return res;
}

ComplexVector ComplexMatrix::GetColumn(int col) const {
    if (col < 0 || col >= cols_)
        throw std::out_of_range("ComplexMatrix::GetColumn - column index out of range");

    ComplexVector v(rows_);

    // data_ 是列主序，列 col 在内存中是连续的： &data_[col * rows_]
    // 可以用 memcpy 提速
    std::memcpy(v.Data(), &data_[static_cast<size_t>(col) * rows_],
                static_cast<size_t>(rows_) * sizeof(MKL_Complex16));

    return v;
}

ComplexMatrix ComplexMatrix::GetRow(const IntArray &row) const {
    ComplexMatrix res(row.giveSize(), cols_);
    for (int i = 0; i < row.giveSize(); ++i) {
        for (int j = 0; j < cols_; ++j) {
            res(i, j) = (*this)(row(i), j);
        }
    }
    return res;
}

void ComplexMatrix::SetSubMatrix(const IntArray &row, const IntArray &col, const ComplexMatrix &sub) {
    for (int i = 0; i < row.giveSize(); ++i) {
        for (int j = 0; j < col.giveSize(); ++j) {
            (*this)(row(i),col(j)) = sub.At(i,j);
        }
    }
}


ComplexMatrix ComplexMatrix::SetColumn(int col, const ComplexVector &value) {
    if (col < 0 || col >= cols_)
        throw std::out_of_range("ComplexMatrix::SetColumn - column index out of range");

    if (value.Size() != rows_)
        throw std::invalid_argument("ComplexMatrix::SetColumn - size mismatch");

    // 列主序：第 col 列从 data_[col * rows_] 开始
    size_t offset = static_cast<size_t>(col) * rows_;

    for (int i = 0; i < rows_; ++i) {
        data_[offset + i] = value[i];   // 直接赋值
    }

    return *this;
}


IntArray ComplexMatrix::FindMaxElementLoaction2() {
    IntArray row(cols_);

    FloatArray Vector1(rows_),Vector2(rows_);
    for (int i=0; i<rows_; ++i) {
        Vector1(i)=this->At(i,0).norm();
    }
    row(0)=Vector1.giveIndexMaxElem()-1;
    auto ele_a=this->At(row(0),1);
    auto ele_b=this->At(row(0),0);
    auto alpha=ele_a/ele_b;
    for (int i=0; i<rows_; ++i) {
        Vector2(i)=(this->At(i,1)-alpha*this->At(i,0)).norm();
    }
    row(1)=Vector2.giveIndexMaxElem()-1;
    return row;
}





// ======================== 初始化 ========================
void ComplexMatrix::Zero()
{
    MKL_Complex16 zero = {0.0, 0.0};
    std::fill(data_.begin(), data_.end(), zero);
}




void ComplexMatrix::findNonZero(IntArray& rows, IntArray& cols) const {
    rows.zero();
    cols.zero();
    for (int j = 0; j < cols_; ++j)
    {
        for (int i = 0; i < rows_; ++i)
        {
            const MKL_Complex16 &v = data_[i + j * rows_];
            if (v.real != 0.0 || v.imag != 0.0)
            {
                rows.append(i);
                cols.append(j);
            }
        }
    }
}

// ======================== 复制构造和赋值 ========================
ComplexMatrix::ComplexMatrix(const ComplexMatrix &other)
    : rows_(other.rows_), cols_(other.cols_), data_(other.data_)
{}

ComplexMatrix &ComplexMatrix::operator=(const ComplexMatrix &other)
{
    if(this != &other)
    {
        rows_ = other.rows_;
        cols_ = other.cols_;
        data_ = other.data_;
    }
    return *this;
}

// ======================== 矩阵运算 ========================
ComplexMatrix ComplexMatrix::operator+(const ComplexMatrix &other) const
{
    if(rows_ != other.rows_ || cols_ != other.cols_)
        throw std::runtime_error("Matrix size mismatch for addition");

    ComplexMatrix res(rows_, cols_);
    for(size_t i = 0; i < data_.size(); i++)
        res.data_[i] = {data_[i].real + other.data_[i].real,
                        data_[i].imag + other.data_[i].imag};
    return res;
}

ComplexMatrix ComplexMatrix::operator-(const ComplexMatrix &other) const
{
    if(rows_ != other.rows_ || cols_ != other.cols_)
        throw std::runtime_error("Matrix size mismatch for subtraction");

    ComplexMatrix res(rows_, cols_);
    for(size_t i = 0; i < data_.size(); i++)
        res.data_[i] = {data_[i].real - other.data_[i].real,
                        data_[i].imag - other.data_[i].imag};
    return res;
}

ComplexMatrix ComplexMatrix::operator*(const ComplexMatrix &other) const
{
    if(cols_ != other.rows_)
        throw std::runtime_error("Matrix size mismatch for multiplication");

    ComplexMatrix res(rows_, other.cols_);

    MKL_Complex16 alpha = {1.0, 0.0};
    MKL_Complex16 beta = {0.0, 0.0};

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            rows_, other.cols_, cols_,
            &alpha,
            data_.data(), rows_,
            other.data_.data(), other.rows_,
            &beta,
            res.data_.data(), res.rows_);

    return res;
}

ComplexMatrix ComplexMatrix::operator*(const ComplexScalar &alpha) const
{
    ComplexMatrix res(*this);
    MKL_Complex16 aval = {alpha.real(), alpha.imag()};
    cblas_zscal(res.Size(), &aval, res.data_.data(), 1);
    return res;
}

// ======================== 矩阵-向量乘法 ========================
ComplexVector ComplexMatrix::operator*(const ComplexVector &vec) const
{
    if(cols_ != vec.Size())
        throw std::runtime_error("Matrix and vector size mismatch");

    ComplexVector res(rows_);
    MKL_Complex16 alpha = {1.0, 0.0};
    MKL_Complex16 beta  = {0.0, 0.0};

    cblas_zgemv(CblasColMajor, CblasNoTrans,
                rows_, cols_,
                &alpha,
                data_.data(), rows_,
                vec.Data(), 1,
                &beta,
                res.Data(), 1);
    return res;
}

// ======================== 输出 ========================
void ComplexMatrix::Print(std::ostream &out) const
{
    int precision = std::numeric_limits<double>::max_digits10;
    for(int i = 0; i < rows_; i++)
    {
        for(int j = 0; j < cols_; j++)
        {
            auto &v = (*this)(i,j);
            out << std::setprecision(precision)
                << v.real
                << " " << (v.imag >= 0 ? "+" : "-")
                << " " << std::abs(v.imag) << "i"<<" ";

        }
        out << "\n";
    }
}

void ComplexMatrix::ComputeEigenValueAndEigenVector(ComplexVector& EigenValue,ComplexMatrix& EigenVector) {
    auto info=LAPACKE_zgeev(CblasColMajor,'N','V',this->Rows(),this->Data(),this->Rows(),EigenValue.Data(),nullptr,this->Rows(),
        EigenVector.Data(),this->Rows());
}

ComplexMatrix ComplexMatrix::ConcatenateRow(const ComplexMatrix& A, const ComplexMatrix& B)
{
    if (A.Rows() != B.Rows())
        throw std::invalid_argument("Row mismatch in ConcatenateRow.");

    int rows = A.Rows();
    int cols = A.Cols() + B.Cols();

    ComplexMatrix C(rows, cols);
    C.Zero();

    // Copy A
    for (int j = 0; j < A.Cols(); ++j)
        for (int i = 0; i < rows; ++i)
            C(i, j) = A(i, j);

    // Copy B
    for (int j = 0; j < B.Cols(); ++j)
        for (int i = 0; i < rows; ++i)
            C(i, j + A.Cols()) = B(i, j);

    return C;
}

ComplexMatrix ComplexMatrix::ConcatenateCol(const ComplexMatrix& A, const ComplexMatrix& B)
{
    if (A.Cols() != B.Cols())
        throw std::invalid_argument("Column mismatch in ConcatenateCol.");

    int rows = A.Rows() + B.Rows();
    int cols = A.Cols();

    ComplexMatrix C(rows, cols);
    C.Zero();

    // Copy A
    for (int j = 0; j < cols; ++j)
        for (int i = 0; i < A.Rows(); ++i)
            C(i, j) = A(i, j);

    // Copy B
    for (int j = 0; j < cols; ++j)
        for (int i = 0; i < B.Rows(); ++i)
            C(i + A.Rows(), j) = B(i, j);

    return C;
}

ComplexMatrix ComplexMatrix::ConcatenateDiag(const ComplexMatrix& A, const ComplexMatrix& B)
{
    int rows = A.Rows() + B.Rows();
    int cols = A.Cols() + B.Cols();

    ComplexMatrix C(rows, cols);
    C.Zero();

    // Copy A
    for (int j = 0; j < A.Cols(); ++j)
        for (int i = 0; i < A.Rows(); ++i)
            C(i, j) = A(i, j);

    // Copy B
    for (int j = 0; j < B.Cols(); ++j)
        for (int i = 0; i < B.Rows(); ++i)
            C(i + A.Rows(), j + A.Cols()) = B(i, j);

    return C;
}
