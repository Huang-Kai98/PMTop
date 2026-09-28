// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#pragma once
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

class DofManager {
public:
  // 用于存储节点到全局自由度的映射
  std::map<int, std::vector<int>> node_to_dofs;

  // 单元到全局自由度的映射
  std::map<int, std::vector<int>> cell_to_global_dofs;

  // 自由度的约束标记（如Dirichlet边界）
  std::map<int, bool> constrained_dofs;

  // CE关联信息
  std::map<int, std::pair<std::vector<int>, std::vector<double>>>
      CE_association;

  // CP关联信息
  std::unordered_map<int, int> CP_association;

  std::vector<int> cp_dofs_list;
  // 悬挂节点约束
  std::map<int, std::vector<int>> hanging_node_constraints;

  // 周期性边界自由度映射
  std::map<int, int> periodic_dofs_mapping;

  DofManager(){};

  // 添加节点自由度映射
  void addNode(int node_id, const std::vector<int> &dofs) {
    node_to_dofs[node_id] = dofs;
  }

  // 添加单元自由度映射
  void addCell(int cell_id, const std::vector<int> &dofs) {
    cell_to_global_dofs[cell_id] = dofs;
  }

  // 设置悬挂节点约束
  void setHangingNodeConstraint(int dof_id,
                                const std::vector<int> &master_dofs) {
    hanging_node_constraints[dof_id] = master_dofs;
  }

  // 设置周期性边界条件
  void setPeriodicBoundary(int dof_master, int dof_slave) {
    periodic_dofs_mapping[dof_slave] = dof_master;
  }

  // 施加Dirichlet边界条件
  void constrainDoF(int dof_id) { constrained_dofs[dof_id] = true; }

  // 施加CE约束
  void addCE(int CE_id, const std::vector<int> &dofs,
             const std::vector<double> &values) {
    CE_association[CE_id] = std::make_pair(dofs, values);
  }

  void addCP(int CP_id, int CE_id) { CP_association[CP_id] = CE_id; }

  // 输出自由度映射
  void printDoFMapping() const {
    std::cout << "Node to DoF mapping:\n";
    for (const auto &node : node_to_dofs) {
      std::cout << "Node " << node.first << ": ";
      for (int dof : node.second) {
        std::cout << dof << " ";
      }
      std::cout << "\n";
    }
  }

  // 处理悬挂节点
  void applyHangingNodeConstraints() {
    for (const auto &constraint : hanging_node_constraints) {
      int slave_dof = constraint.first;
      const std::vector<int> &master_dofs = constraint.second;
      // 这里可以定义如何将悬挂节点的自由度与主节点自由度联系起来
      // 例如，计算加权平均
      std::cout << "Applying hanging node constraint for DoF " << slave_dof
                << std::endl;
    }
  }

  // 处理周期性边界条件
  void applyPeriodicConstraints() {
    for (const auto &periodic_pair : periodic_dofs_mapping) {
      int slave_dof = periodic_pair.first;
      int master_dof = periodic_pair.second;
      std::cout << "DoF " << slave_dof << " is periodic with DoF " << master_dof
                << std::endl;
    }
  }

  // 获取总体自由度数
  int getTotalDoFs() const {
    std::set<int> unique_dofs;

    // 遍历所有节点，收集自由度
    for (const auto &node : node_to_dofs) {
      unique_dofs.insert(node.second.begin(), node.second.end());
    }

    // 遍历所有单元，收集自由度
    for (const auto &cell : cell_to_global_dofs) {
      unique_dofs.insert(cell.second.begin(), cell.second.end());
    }

    return unique_dofs.size(); // 返回去重后的自由度总数
  }
};
