#pragma once
#include <iostream>
#include <ostream>
#include <cmath>

class Node {
public:
  int id;
  double x, y, z, thxy, thyz, thzx;
  Node(int id_, double x_, double y_, double z_,double thxy_,double thyz_,double thzx_):
  id(id_),x(x_),y(y_),z(z_),thxy(thxy_),thyz(thyz_),thzx(thzx_){}

  Node() : id(0), x(0), y(0), z(0), thxy(0), thyz(0), thzx(0) {}

  Node(const Node &other)
      : id(other.id), x(other.x), y(other.y), z(other.z), thxy(other.thxy),
        thyz(other.thyz), thzx(other.thzx) {}
  // 定义友元函数用于输出
  friend std::ostream &operator<<(std::ostream &os, const Node &node) {
    os << "(" << node.id << ", " << node.x << ", " << node.y << ", " << node.z
       << ", " << node.thxy << ", " << node.thyz << ", " << node.thzx << ")"
       << std::endl;
    return os;
  }

  double euclideanDistance(const Node& other) {
    double dx = other.x - x;
    double dy = other.y - y;
    double dz = other.z - z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  void printNode() {
    std::cout << "id:" << id << " " << x << " " << y << " " << z << " " << thxy
              << " " << thyz << " " << thzx << std::endl;
  }
};