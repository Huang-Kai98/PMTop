// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include <mkl.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <regex>
#include <fstream>
#include <iomanip>

#include "ComplexMatrix.h"
#include "ComplexScalar.hpp"
#include "floatarray.h"
#include "intarray.h"

class ComplexMatrix;

/// 使用 MKL COMPLEX16 内部函数实现的复数向量类
class ComplexVector {
private:
    std::vector<MKL_Complex16> data;

public:
    /// 默认构造
    ComplexVector() = default;

    /// 构造长度为 n 的向量并初始化为 0
    explicit ComplexVector(int n) : data(n) {
        for (auto &z : data) { z.real = 0.0; z.imag = 0.0; }
    }

    /// 从 std::vector<std::complex<double>> 构造
    explicit ComplexVector(const std::vector<std::complex<double>> &v) {
        data.resize(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            data[i].real = v[i].real();
            data[i].imag = v[i].imag();
        }
    }

    ComplexVector(const FloatArray &real_part,
                     const FloatArray &imag_part) {
        assert(real_part.giveSize() == imag_part.giveSize());
        size_t n = real_part.giveSize();
        data.resize(n);
        for (size_t i = 0; i < n; ++i) {
            data[i].real = real_part[i];
            data[i].imag = imag_part[i];
        }
    }

    void EliminateElement(const IntArray &index) {
        for (int i=0; i<index.giveSize(); ++i) {
                this->Set(index(i),0.0,0.0);
        }
    }
    /// 拷贝构造
    ComplexVector(const ComplexVector &other) = default;

    ComplexVector(std::string file) {
        std::ifstream fin(file);
        if (!fin) {   // ← 添加文件打开失败检查
            std::cerr << "Error: failed to open eigenvector.txt" << std::endl;
        }
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty()) continue;
            this->push_back(parseComplex(line));
        }
        fin.close();
    }

    void push_back(const ComplexScalar &val) {
        data.push_back(val);
    }


    /// 赋值
    ComplexVector &operator=(const ComplexVector &other) = default;

    /// 返回向量长度
    int Size() const { return static_cast<int>(data.size()); }

    /// 元素访问
    ComplexScalar operator()(int i) const {
        assert(i >= 0 && i < Size());
        return data[i];
    }

    IntArray findNonZero() const {
        IntArray res;
        res.resize(Size());
        for (int i = 0; i < Size(); ++i)
        {
            const MKL_Complex16 &v = data[i];
            if (v.real != 0.0 || v.imag != 0.0)
                res.append(i);
        }
        return res;
    }

    // 解析形如 "a + bi" 或 "a - bi" 的复数字符串
    ComplexScalar parseComplex(const std::string& s) {
        const char* p = s.c_str();

        // 1. 解析实部
        char* endptr;
        double real = std::strtod(p, &endptr);

        // 移动到符号位置 ('+' 或 '-')
        p = endptr;
        while (*p == ' ') p++;

        char sign = *p; // '+' 或 '-'
        p++;

        // 跳过空格
        while (*p == ' ') p++;

        // 2. 解析虚部（后面会有一个 'i'）
        double imag = std::strtod(p, &endptr);

        // 根据符号调整虚部符号
        if (sign == '-') imag = -imag;

        return {real, imag};

    }

    /// 只读访问
    const ComplexScalar operator[](int i) const {
        assert(i >= 0 && i < Size());
        return data[i];
    }


    MKL_Complex16 &operator()(int i)
    {
        return data[i];
    }

    ComplexVector operator()(IntArray &index) {
        ComplexVector res(index.giveSize());
        for (int i=0;i<index.giveSize();i++) {
            res.Set(i,(*this)(index[i]));
        }
        return res;
    }

    ComplexVector GetRow(IntArray &index) {
        return (*this)(index);
    }

    ComplexVector Conjugate() {
        ComplexVector res(*this);
        for (int i = 0; i < Size(); i++) {
            res.Set(i,data[i].real,-data[i].imag);
        }
        return res;
    }


    int FindMinElementLoaction() {
        FloatArray Vector1(Size());
        for (int i=0; i<Size(); ++i) {
            Vector1(i)=sqrt(this->data[i].real*this->data[i].real+this->data[i].imag*this->data[i].imag);
        }


        return Vector1.giveIndexMinElem()-1;;
    }

    int FindMaxElementLoaction() {
        FloatArray Vector1(Size());
        for (int i=0; i<Size(); ++i) {
            Vector1(i)=sqrt(this->data[i].real*this->data[i].real+this->data[i].imag*this->data[i].imag);
        }


        return Vector1.giveIndexMaxElem()-1;;
    }


    /// 使用 ComplexScalar 设置元素
    void Set(int i, const ComplexScalar &val) {
        assert(i >= 0 && i < Size());
        data[i].real = val.real();
        data[i].imag = val.imag();
    }

    /// 使用实部和虚部设置元素
    void Set(int i, double real_part, double imag_part) {
        assert(i >= 0 && i < Size());
        data[i].real = real_part;
        data[i].imag = imag_part;
    }

    /// 全部置零
    void Zero() {
        MKL_Complex16 zero = {0.0, 0.0};
        cblas_zscal(Size(), &zero, data.data(), 1);
    }

    /// 向量加法: this + rhs
    ComplexVector operator+(const ComplexVector &rhs) const {
        assert(Size() == rhs.Size());
        ComplexVector res(*this);
        MKL_Complex16 one = {1.0, 0.0};
        cblas_zaxpy(Size(), &one, rhs.data.data(), 1, res.data.data(), 1);
        return res;
    }

    /// 向量减法: this - rhs
    ComplexVector operator-(const ComplexVector &rhs) const {
        assert(Size() == rhs.Size());
        ComplexVector res(*this);
        MKL_Complex16 minus_one = {-1.0, 0.0};
        cblas_zaxpy(Size(), &minus_one, rhs.data.data(), 1, res.data.data(), 1);
        return res;
    }

    /// 向量标量乘法
    ComplexVector operator*(const ComplexScalar &alpha) const {
        ComplexVector res(*this);
        MKL_Complex16 a = {alpha.real(), alpha.imag()};
        cblas_zscal(Size(), &a, res.data.data(), 1);
        return res;
    }

    /// 向量标量乘法
    ComplexVector operator*(const double &alpha) const {
        ComplexVector res(*this);
        MKL_Complex16 a = {alpha, 0.0};
        cblas_zscal(Size(), &a, res.data.data(), 1);
        return res;
    }

    /// 向量向量乘法
    ComplexScalar operator*(const ComplexVector &alpha) const {
        return this->Dot(alpha);  // 返回复数标量
    }

    /// 向量标量除法
    ComplexVector operator/(const ComplexScalar &alpha) const {
        ComplexVector res(*this);
        std::complex<double> inv = 1.0 / static_cast<std::complex<double>>(alpha);
        MKL_Complex16 a = {inv.real(), inv.imag()};
        cblas_zscal(Size(), &a, res.data.data(), 1);
        return res;
    }

    /// 复共轭内积 dot = conj(this)^T * rhs
    ComplexScalar conjDot(const ComplexVector &rhs) const {
        assert(Size() == rhs.Size());
        MKL_Complex16 result;
        cblas_zdotc_sub(Size(), data.data(), 1, rhs.data.data(), 1, &result);
        return ComplexScalar(result.real, result.imag);
    }

    /// 复内积 dot = (this)^T * rhs
    ComplexScalar Dot(const ComplexVector &rhs) const {
        assert(Size() == rhs.Size());
        MKL_Complex16 result;
        cblas_zdotu_sub(Size(), data.data(), 1, rhs.data.data(), 1, &result);
        return ComplexScalar(result.real, result.imag);
    }

    /// 二范数 ||x||
    double Norm2() const {
        return cblas_dznrm2(Size(), data.data(), 1);
    }

    /// 判断是否为空
    bool empty() const {
        return data.empty();
    }

    /// 调整向量长度并初始化（默认全零，可指定初值）
    void resize(int n, const ComplexScalar &init_val = ComplexScalar(0.0, 0.0)) {
        data.resize(n);
        MKL_Complex16 val;
        val.real = init_val.real();
        val.imag = init_val.imag();
        std::fill(data.begin(), data.end(), val);
    }

    /// 打印
    void Print(std::ostream &out = std::cout) const {
        int precision = std::numeric_limits<double>::max_digits10;
        for (int i = 0; i < Size(); ++i) {
            auto z = data[i];
            out << std::setprecision(precision)
                << z.real
                << " " << (z.imag >= 0 ? "+" : "-")
                << " " << std::abs(z.imag) << "i";

            if (i < Size() - 1)
                out << "\n";
        }
        out << "\n";
    }


    /// 返回底层指针（可直接用于 MKL 接口）
    MKL_Complex16 *Data() { return data.data(); }
    const MKL_Complex16 *Data() const { return data.data(); }
};


inline std::ostream &operator<<(std::ostream &out, const ComplexVector &v)
{
    for (int i = 0; i < v.Size(); ++i)
    {
        auto z = v(i);

        out << z.real()
            << " " << (z.imag() >= 0 ? "+" : "-")
            << " " << std::abs(z.imag()) << "i";

        if (i < v.Size() - 1)
            out << "\n";
    }
    return out;
}