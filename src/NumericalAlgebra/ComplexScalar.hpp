// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include <mkl.h>
#include <complex>
#include <iostream>

/// 封装 MKL COMPLEX16 的复数标量类
class ComplexScalar {
private:
    MKL_Complex16 value;

public:
    /// 默认构造：0 + 0i
    ComplexScalar() {
        value.real = 0.0;
        value.imag = 0.0;
    }

    ComplexScalar(const MKL_Complex16& z) {
        value.real = z.real;
        value.imag = z.imag;
    }

    MKL_Complex16 givePointer() const {
        return value;
    }


    // ✅ 添加转换运算符，支持直接赋值给 MKL_Complex16
    operator MKL_Complex16() const {
        MKL_Complex16 z;
        z.real = value.real;
        z.imag = value.imag;
        return z;
    }

    /// 构造：实部和虚部
    ComplexScalar(double re, double im = 0.0) {
        value.real = re;
        value.imag = im;
    }

    /// 从 std::complex<double> 构造
    ComplexScalar(const std::complex<double> &z) {
        value.real = z.real();
        value.imag = z.imag();
    }

    /// 拷贝构造
    ComplexScalar(const ComplexScalar &other) {
        value = other.value;
    }

    /// 赋值
    ComplexScalar &operator=(const ComplexScalar &other) {
        value = other.value;
        return *this;
    }

    /// 转换为 std::complex<double>
    operator std::complex<double>() const {
        return std::complex<double>(value.real, value.imag);
    }

    /// 获取实部
    double real() const { return value.real; }
    double& real() { return value.real; }       // 非 const 版本，可赋值

    /// 获取虚部
    double imag() const { return value.imag; }
    double& imag() { return value.imag; }       // 非 const 版本，可赋值

    /// 设置数值
    void set(double re, double im = 0.0) {
        value.real = re;
        value.imag = im;
    }

    /// 基本运算符
    ComplexScalar operator+(const ComplexScalar &rhs) const {
        return ComplexScalar(value.real + rhs.value.real,
                             value.imag + rhs.value.imag);
    }

    ComplexScalar operator-(const ComplexScalar &rhs) const {
        return ComplexScalar(value.real - rhs.value.real,
                             value.imag - rhs.value.imag);
    }

    ComplexScalar operator*(const ComplexScalar &rhs) const {
        double re = value.real * rhs.value.real - value.imag * rhs.value.imag;
        double im = value.real * rhs.value.imag + value.imag * rhs.value.real;
        return ComplexScalar(re, im);
    }

    ComplexScalar operator/(const ComplexScalar &rhs) const {
        double denom = rhs.value.real * rhs.value.real + rhs.value.imag * rhs.value.imag;
        double re = (value.real * rhs.value.real + value.imag * rhs.value.imag) / denom;
        double im = (value.imag * rhs.value.real - value.real * rhs.value.imag) / denom;
        return ComplexScalar(re, im);
    }

    ComplexScalar &operator+=(const ComplexScalar &rhs) {
        value.real += rhs.value.real;
        value.imag += rhs.value.imag;
        return *this;
    }

    ComplexScalar &operator-=(const ComplexScalar &rhs) {
        value.real -= rhs.value.real;
        value.imag -= rhs.value.imag;
        return *this;
    }

    ComplexScalar &operator*=(const ComplexScalar &rhs) {
        double re = value.real * rhs.value.real - value.imag * rhs.value.imag;
        double im = value.real * rhs.value.imag + value.imag * rhs.value.real;
        value.real = re;
        value.imag = im;
        return *this;
    }

    ComplexScalar &operator/=(const ComplexScalar &rhs) {
        double denom = rhs.value.real * rhs.value.real + rhs.value.imag * rhs.value.imag;
        double re = (value.real * rhs.value.real + value.imag * rhs.value.imag) / denom;
        double im = (value.imag * rhs.value.real - value.real * rhs.value.imag) / denom;
        value.real = re;
        value.imag = im;
        return *this;
    }

    ComplexScalar sqrt() {
        ComplexScalar res;
        vzSqrt(1, &value, &res.value);  // 计算复数开根号
        return res;
    }

    double norm() const {
        return std::sqrt(value.real * value.real + value.imag * value.imag);
    }

    /// 打印
    friend std::ostream &operator<<(std::ostream &os, const ComplexScalar &z) {
        os << "(" << z.value.real << (z.value.imag >= 0 ? "+" : "") << z.value.imag << "i)";
        return os;
    }
};
