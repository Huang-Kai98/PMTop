#pragma once
#include "mem_alloc.hpp"

// 存储三维关系
// STable3D 用于存储三维网格或实体之间的稀疏连接关系。
// 每个连接由三维索引 (r, c, f)（行、列、层）表示，可能额外扩展到四维 (r, c, f,
// t)。
// 表示的三维关系是对称的，例如 (r, c, f) 等价于 (c, r, f)。
// 对称性允许减少存储和计算成本。
// 只存储非零连接，使用链表结构 STable3DNode 实现。
// 支持高效地插入和查询连接。
//

class STable3DNode {
public:
  STable3DNode *Prev;
  int Column, Floor, Number;
};

/// Symmetric 3D Table
class STable3D {
private:
  int Size, NElem;
  STable3DNode **Rows;

#ifdef MFEM_USE_MEMALLOC
  MemAlloc<STable3DNode, 1024> NodesMem;
#endif

public:
  STable3D(int nr);

  int Push(int r, int c, int f);

  int operator()(int r, int c, int f) const;

  int Index(int r, int c, int f) const;

  int Push4(int r, int c, int f, int t);

  int operator()(int r, int c, int f, int t) const;

  int NumberOfElements() { return NElem; };

  ~STable3D();
};