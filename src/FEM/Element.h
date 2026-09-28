#pragma once
#include "Node.h"
#include <iostream>
#include <ostream>
#include <vector>

class Element {
public:
  double *elementStif;
  double *elementMass;
  double *elementDamp;
  int dof;
  int nr;
  std::string jtyp;
  std::vector<Node> nodes;
  int id;
  std::vector<int> MatrixLocaltion;
  bool isTop;
  std::vector<int> keypoint;

  virtual ~Element(){}

  virtual void init() = 0;
  virtual void setStiff(const int i, const int j, const double value) = 0;
  virtual void setStiff(const int i, const double value) {
    elementStif[i] = value;
  }
  virtual void setMass(const int i, const int j, const double value) = 0;
  virtual void setMass(const int i, const double value) {
    elementMass[i] = value;
  }

  virtual void setDamp(const int i, const int j, const double value) = 0;
  virtual void setDamp(const int i, const double value) {
    elementDamp[i] = value;
  }


  double getStiff(const int i, const int j) const {
    return elementStif[j * nr + i];
  }

  double getMass(const int i, const int j) const {
    return elementMass[j * nr + i];
  }

  double getDamp(const int i, const int j) const {
    return elementDamp[j * nr + i];
  }

  void setMatrixLoacltion() {

    MatrixLocaltion.clear();
    for (auto const &node : nodes) {
      if (jtyp == "COMBIN14") {
        if (keypoint[1]==0) {
          MatrixLocaltion.push_back((node.id-1)*6+1);
          MatrixLocaltion.push_back((node.id-1)*6+2);
          MatrixLocaltion.push_back((node.id-1)*6+3);
        }
        else {
          MatrixLocaltion.push_back((node.id-1)*6+keypoint[1]);
        }
      }
      else {
        for (int j = 0; j < dof; ++j) {
          MatrixLocaltion.push_back((node.id-1)*dof+j+1);
        }
      }

    }


    // int size=AnasysElements[i]->nr;
    // std::vector<Node> nodes = AnasysElements[i]->nodes;
    // AnasysElements[i]->MatrixLocaltion.clear();
    // for (auto const &node : nodes) {
    //   std::vector<int> tempnode=dof_manager.node_to_dofs[node.id];
    //   for (auto const &ee:tempnode) {
    //     AnasysElements[i]->MatrixLocaltion.push_back(ee);
    //   }
    // }
  }
  int getDof() const { return dof; }

  int getNr() const { return nr; }

  void setKeypoint(const int* keyopts) {
    keypoint.resize(18);
    for (int i = 0; i < 18; i++) {
      keypoint[i] = keyopts[i];
    }
  }

  void setDof(const int value) { dof = value; }

  void setNr(const int value) { nr = value; }

  void setJtyp(const std::string type) { jtyp = type; }

  void setNodes(const Node node) { nodes.push_back(node); }
  void saveToStream(std::ostream &out) const {
    // 输出 elementStif
    out << "Element Stiffness: " << std::endl;
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nr; j++) {
        out << elementStif[j * nr + i] << " ";
      }
      out << std::endl;
    }
    out << std::endl;

    // 输出 elementMass
    out << "Element Mass: " << std::endl;
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nr; j++) {
        out << elementMass[j * nr + i] << " ";
      }
      out << std::endl;
    }
    out << std::endl;

    // 输出 elementMass
    out << "Element Damp: " << std::endl;
    for (int i = 0; i < nr; ++i) {
      for (int j = 0; j < nr; j++) {
        out << elementDamp[j * nr + i] << " ";
      }
      out << std::endl;
    }
    out << std::endl;

    // 输出 dof, nr, jtyp, id
    out << "DOF: " << dof << std::endl;
    out << "NR: " << nr << std::endl;
    out << "JType: " << jtyp << std::endl;
    out << "ID: " << id << std::endl;

    // 输出 nodes
    out << "Nodes: " << std::endl;
    for (const Node &node : nodes) {
      out << node;
    }
    out << std::endl;

    // 输出 MatrixLocaltion
    out << "Matrix Locations: " << std::endl;
    for (int loc : MatrixLocaltion) {
      out << loc << " ";
    }
    out << "                      " << std::endl;
    out << std::endl;
  }

  void setID(const int value) { id = value; }
  void printElement() {
    std::cout << "dof:" << dof << std::endl;
    std::cout << "nr:" << nr << std::endl;
    std::cout << "jtyp:" << jtyp << std::endl;
    std::cout << "id:" << id << std::endl;
    for (auto node : nodes) {
      std::cout << "node:" << node.id << std::endl;
    }
    std::cout << "elementStif:" << std::endl;
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nr; j++) {
        std::cout << elementStif[j * nr + i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "elementMass:" << std::endl;
    for (int i = 0; i < nr; i++) {
      for (int j = 0; j < nr; j++) {
        std::cout << elementMass[j * nr + i] << " ";
      }
      std::cout << std::endl;
    }
  }

  void setTop(const bool status){ isTop = status; };
  bool getTopStatus() const{return isTop;};

};