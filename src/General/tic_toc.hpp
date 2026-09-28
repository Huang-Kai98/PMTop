// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once

#include <sys/times.h>

/// Timing object
// 成员变量：
// real_time, user_time, syst_time:
// 分别记录实际时间、用户时间和系统时间的累积值。
// start_rtime, start_utime, start_stime:
// 分别记录开始计时时的实际时间、用户时间和系统时间。
// my_CLK_TCK:
// 表示每秒钟的时钟计时单位（如每秒 CLOCKS_PER_SEC 个时钟周期）。
// Running:
// 标记计时器当前是否正在运行，避免多次启动和停止时计算错误.
// Current:
// 私有方法，用于获取当前时间信息（需要结合 times() 或 clock() 函数实现）。
// StopWatch():
// 构造函数，初始化时钟频率 my_CLK_TCK.
// Clear():
// 重置所有计时值为零。
// Start() 和 Stop():
// 分别开始和停止计时，更新时间差。
// RealTime(), UserTime(), SystTime():
// 分别返回累计的实际时间、用户时间和系统时间。
// 全局辅助函数
// tic() 和 toc():
// 通过使用全局 StopWatch 对象（tic_toc）实现快速计时的启动和获取间隔时间。
class StopWatch {
private:
  clock_t real_time, user_time, syst_time;
  clock_t start_rtime, start_utime, start_stime;
  long my_CLK_TCK;
  short Running;
  void Current(clock_t *, clock_t *, clock_t *);

public:
  StopWatch(); // determines my_CLK_TCK
  void Clear();
  void Start();
  void Stop();
  double RealTime();
  double UserTime();
  double SystTime();
};

extern StopWatch tic_toc;

/// Start timing
extern void tic();

/// End timing
extern double toc();