#pragma once
#include "array.hpp"
#include "table.hpp"

// IntegerSet 提供了基础的整数集合操作，适合用于管理简单的整数分组。
// ListOfIntegerSets 管理多个整数集合的操作，例如插入、查找等。
// 借助 AsTable，可以将整数集合的关系表示为稀疏表（Table
// 类型），用于科学计算或图结构操作。
//

// 类 IntegerSet
// 表示一个整数集合，封装了一些常用的集合操作。
// Size()：返回集合的大小。
// operator Array<int>&()：支持将 IntegerSet 转换为 Array<int> 类型
// PickElement()：选择集合中的第一个元素。
// PickRandomElement()：随机选择集合中的一个元素。
// operator==(IntegerSet &s)：判断两个集合是否相等。
// Recreate(const int n, const int *p)：用新的元素重新构建集合。
//
// me：一个 Array<int> 类型对象，存储集合的整数
//

/// A set of integers
class IntegerSet {
private:
  Array<int> me;

public:
  IntegerSet() {}

  IntegerSet(IntegerSet &s);

  IntegerSet(const int n, const int *p) { Recreate(n, p); }

  int Size() { return me.Size(); }

  operator Array<int> &() { return me; }

  int PickElement() { return me[0]; }

  int PickRandomElement();

  int operator==(IntegerSet &s);

  void Recreate(const int n, const int *p);
};
// 类 ListOfIntegerSets
// 用于管理多个 IntegerSet 的列表。
// 提供集合插入、查找等功能。
//
// Size()：返回列表中集合的数量。
// PickElementInSet(int i)：从第i个集合中选择第一个元素。
// PickRandomElementInSet(int i)：从第i个集合中随机选择一个元素。
// Insert(IntegerSet &s)：插入一个新的整数集合。
// Lookup(IntegerSet &s)：查找一个整数集合在列表中的位置。
// AsTable(Table &t)：将列表的集合信息转化为 Table 表的形式。
//
// TheList：一个 Array<IntegerSet *> 对象，动态存储整数集合的指针。

/// List of integer sets
class ListOfIntegerSets {
private:
  Array<IntegerSet *> TheList;

public:
  int Size() { return TheList.Size(); }

  int PickElementInSet(int i) { return TheList[i]->PickElement(); }

  int PickRandomElementInSet(int i) { return TheList[i]->PickRandomElement(); }

  int Insert(IntegerSet &s);

  int Lookup(IntegerSet &s);

  void AsTable(Table &t);

  ~ListOfIntegerSets();
};