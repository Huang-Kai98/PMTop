// PMTop - finite element and modal sensitivity research code.
// Project maintainer: Huang Kai.
// Copyright (c) 2026 Huang Kai, for original PMTop contributions.
// PMTop contributions are licensed under LGPL-2.1-only; see LICENSE.
// Existing third-party notices and terms remain applicable; see NOTICE.md.

#include "AnsysLoad.h"

#include <algorithm>
#include <cmath>
#include <floatarray.h>

#include "Shell.h"
#include "Solid.h"
#include "Beam.h"
#include "vtk.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
void AnsysLoad::ReadAnsysElementType() {
  AnasysElementType.readFile(filename + "_ETLIST.lis");
  LOG_INFO("ElementTypeReader Success!");
  LOG_INFO("Number of ElementTypes: {}",AnasysElementType.getElementCount());
}

void AnsysLoad::ReadAnsysNode() {
  std::ifstream file(filename + "_NLIST.lis");
  std::string line;
  while (std::getline(file, line)) {
    if (line.find("NODE") != std::string::npos) {
      continue; // Skip header lines
    }
    std::istringstream iss(line);
    Node node;
    if (iss >> node.id >> node.x >> node.y >> node.z >> node.thxy >>
        node.thyz >> node.thzx) {
      AnasysNodes.push_back(node);
    }
  }
  file.close();
  LOG_INFO("Node Reader Success!");
  LOG_INFO("Number of Nodes: {}", AnasysNodes.size());
  NumberOfNode = AnasysNodes.size();
}

void AnsysLoad::ReadAnasysElementNodes() {
  std::ifstream file(filename + "_ELIST.lis");
  if (!file.is_open()) {
    LOG_ERROR("Error opening file!");
  }
  std::string line;
  // 跳过文件的标题行
  while (std::getline(file, line)) {
    if (line.find("ELEM") != std::string::npos) {
      break; // 找到标题行后停止跳过
    }
  }
  // 读取元素数据
  while (std::getline(file, line)) {
    if (line.empty())
      continue; // 跳过空行
    ElementNodesInfo element = parseLine(line);

    // 只保存那些有有效元素编号且节点数量大于0的元素
    if (element.elem > 0 && !element.nodes.empty()) {
      elements[element.elem] = element; // 使用elem作为map的键
    }
  }
  file.close();
  NumberOfElement = elements.size();
  LOG_INFO("ElementNodesInfo Success!");
  LOG_INFO("Number of Elements: {}",elements.size());
}

ElementNodesInfo AnsysLoad::parseLine(const std::string &line) {
  std::istringstream iss(line);
  ElementNodesInfo data;
  int node;

  // 初始化元素编号等为0
  data.elem = data.mat = data.typ = data.rel = data.esy = data.sec = 0;

  // 读取前6列：ELEM, MAT, TYP, REL, ESY, SEC
  if (!(iss >> data.elem >> data.mat >> data.typ >> data.rel >> data.esy >>
        data.sec)) {
    // 如果未能成功读取前6列，则返回无效的 ElementData
    return data;
  }

  // 读取NODES列，假设每行有任意数量的节点
  while (iss >> node) {
    data.nodes.push_back(node);
  }

  return data;
}

void AnsysLoad::ReadAnasysBoundaryCondition() {
  std::ifstream inputFile(filename + "_DLIST.lis"); // 替换为你的文件名
  std::string line;

  if (!inputFile.is_open()) {
    LOG_ERROR("无法打开文件DLIST!");
    throw std::runtime_error("无法打开文件!");
  }

  // 读取标题行
  while (std::getline(inputFile, line)) {
    if (line.find("CURRENTLY SELECTED DOF SET") != std::string::npos) {
      // 解析 DOF 标签
      std::istringstream iss(line);
      std::string temp;
      // 跳过前面的部分
      while (iss >> temp) {
        if (temp == "CURRENTLY")
          continue;
        if (temp == "SELECTED")
          continue;
        if (temp == "DOF")
          continue;
        if (temp == "SET=")
          continue;
        BoundaryAnasysdofLabels.push_back(temp); // 存储 DOF 标签
      }
      break; // 读取到 DOF 设置后退出循环
    }
  }

  // 读取剩余数据
  while (std::getline(inputFile, line)) {
    // 过滤掉空行和包含 "NODE" 的行
    if (line.empty() || line.find("NODE") != std::string::npos)
      continue;

    std::istringstream iss(line);
    BoundaryData boundaryData;

    // 读取数据
    if (iss >> boundaryData.nodeId >> boundaryData.dofLabel >>
        boundaryData.realPart >> boundaryData.imagPart) {
      BoundaryAnasysData.push_back(boundaryData);
    }
  }

  inputFile.close();

  LOG_INFO("BoundaryReader Success!");
}

void AnsysLoad::ReadCEelement() {
  CE_enable = true;
  std::ifstream file(filename + "_CELIST.lis");
  if (!file.is_open()) {
    LOG_INFO("NO CELIST.");
    CE_enable = false;
  }

  if (CE_enable == true) {
    std::string line;

    while (std::getline(file, line)) {
      if (line.find("CONSTRAINT EQUATION NO.") != std::string::npos) {
        ConstraintEquation eq;
        sscanf(line.c_str(),
               " CONSTRAINT EQUATION NO. %d HAS %d TERMS. CONSTANT= %lf",
               &eq.equation_number, &eq.term_count, &eq.constant);

        for (int i = 0; i < eq.term_count; ++i) {
          Term term;
          char direction[10]; // 临时字符数组用于读取方向

          std::getline(file, line);
          sscanf(line.c_str(), " NODE= %d DIR= %s COEFFICIENT= %lf", &term.node,
                 direction, &term.coefficient);

          term.direction = direction; // 将方向赋值给 std::string 成员
          eq.terms.push_back(term);
        }

        equations.push_back(eq);
      }
    }
    file.close();
    LOG_TRACE("Coupled Sets Read Success!");
    LOG_TRACE("Number of Coupled Sets: {}",equations.size());

  }
}

void AnsysLoad::ReadCoupledSets() {
  CP_enable = true;
  std::ifstream file(filename + "_CPLIST.lis");
  if (!file.is_open()) {
    LOG_INFO("NO CPLIST.");
    CP_enable = false;
  }
  std::string line;

  if (CP_enable) {

    while (std::getline(file, line)) {
      if (line.find("COUPLED SET=") != std::string::npos) {
        CoupledSet set;
        std::istringstream iss(line);
        std::string temp;

        // Parse fields by reading through the line
        iss >> temp >> temp >>
            set.set_id;               // Skip "COUPLED SET=" and read set_id
        iss >> temp >> set.direction; // Skip "DIRECTION=" and read direction
        iss >> temp >> temp >>
            set.total_nodes; // Skip "TOTAL NODES=" and read total_nodes

        // Read the next line for node list
        std::getline(file, line);
        std::istringstream nodeStream(line);
        while (nodeStream >> temp) {   // Skip "NODES=" and read nodes
          if (std::isdigit(temp[0])) { // Ensure we're only reading numbers
            set.nodes.push_back(std::stoi(temp));
          }
        }

        coupledSets.push_back(set);
      }
    }
    LOG_TRACE("Coupled Sets Read Success!");
    LOG_TRACE("Number of Coupled Sets: {}",coupledSets.size());
  }
}

void AnsysLoad::Adjacency_list_node() {
  for (const auto &element : elements) {
    const auto &nodes = element.second.nodes;
    for (size_t i = 0; i < nodes.size(); ++i) {
      for (size_t j = i + 1; j < nodes.size(); ++j) {
        adjacency_list[nodes[i]].insert(nodes[j]);
        adjacency_list[nodes[j]].insert(nodes[i]);
      }
    }
  }
  // 统计顶点和边的总数
  vertex_count = adjacency_list.size();
  edge_count = 0;

  for (const auto &pair : adjacency_list) {
    edge_count += pair.second.size();
  }
  edge_count /= 2; // 因为无向图中的每条边被存储了两次

}

void AnsysLoad::Adjacency_list_element() {
  std::map<std::pair<int, int>, std::set<int>> edge_to_cells; // 边到单元的映射
  for (const auto &cell : elements) {
    const auto &vertices = cell.second.nodes;
    int num_vertices = vertices.size();

    for (int i = 0; i < num_vertices; ++i) {
      for (int j = i + 1; j < num_vertices; ++j) {
        int v1 = std::min(vertices[i], vertices[j]);
        int v2 = std::max(vertices[i], vertices[j]);

        std::pair<int, int> edge = {v1, v2};
        edge_to_cells[edge].insert(cell.second.elem);
      }
    }
  }

  // 根据边和单元的映射，找到所有共享边的单元并记录为邻接单元
  for (const auto &[edge, cell_set] : edge_to_cells) {
    if (cell_set.size() > 1) {
      for (int cell_id : cell_set) {
        for (int neighbor_id : cell_set) {
          if (cell_id != neighbor_id) {
            adjacent_cells[cell_id].insert(neighbor_id);
          }
        }
      }
    }
  }
}

void AnsysLoad::Adjacency_list_element_node() {
  std::map<int, std::set<int>> node_to_cells; // 节点到单元的映射

  // 构建节点到单元的映射
  for (const auto &cell : elements) {
    const auto &vertices = cell.second.nodes;

    for (int vertex : vertices) {
      node_to_cells[vertex].insert(cell.second.elem);
    }
  }

  // 根据节点和单元的映射，找到所有共享节点的单元并记录为邻接单元
  for (const auto &[node, cell_set] : node_to_cells) {
    for (int cell_id : cell_set) {
      for (int neighbor_id : cell_set) {
        if (cell_id != neighbor_id) {
          adjacent_cells[cell_id].insert(neighbor_id);
        }
      }
    }
  }
}

std::unique_ptr<Element> AnsysLoad::createElement(const std::string &type) {
  if (type == "Shell181" || type == "SHELL181") {
    return std::make_unique<Shell181>();
  } else if (type == "Shell63" || type == "SHELL63") {
    return std::make_unique<Shell63>();
  } else if (type == "Beam188" || type == "BEAM188") {
    return std::make_unique<Beam188>();
  }else if (type =="Solid185"|| type == "SOLID185") {
    return std::make_unique<Solid185>();
  }else if (type == "COMBIN14"){
    return std::make_unique<COMBIN14>();
  } else {
    throw std::runtime_error("Unknown element type");
  }
  return nullptr;
}


void AnsysLoad::dofManage() {
  std::map<int, std::set<int>> NodeType;
  std::map<int, int> newNodeType;
  std::map<int, int> elementType;
  for (const auto &element : elements) {
    AnasysElements[element.first]=createElement(AnasysElementType.getElement(element.second.typ)->elementName);
    AnasysElements[element.first]->init();
    AnasysElements[element.first]->setID(element.first);
    AnasysElements[element.first]->setKeypoint(AnasysElementType.getElement(element.second.typ)->keyopts);
    elementType[element.second.typ] = AnasysElements[element.first]->dof;
    for (const auto &node : element.second.nodes) {
      if (node == 0)
        continue;
      AnasysElements[element.first]->setNodes(AnasysNodes[node - 1]);
      NodeType[element.second.typ].insert(node);
      nodeDof[node] = elementType[element.second.typ];
    }
  }






  for (const auto &element : elements) {
    unsigned int nodedof = 0;

    if (AnasysElementType.getElement(element.second.typ)->elementName=="COMBIN14") {
      switch (AnasysElementType.getElement(element.second.typ)->keyopts[1]) {
        case 0: nodedof|=UX|UY|UZ; break;
        case 1: nodedof|=UX; break;
        case 2: nodedof|=UY; break;
        case 3: nodedof|=UZ; break;
        case 4: nodedof|=ROTX; break;
        case 5: nodedof|=ROTY; break;
        case 6: nodedof|=ROTZ; break;
        default: std::cerr << "Unknown element type " << std::endl; break;
      }
    }
    else {

      switch (AnasysElements[element.first]->dof) {
        case 6: nodedof|=UX|UY|UZ|ROTX|ROTY|ROTZ; break;
        case 3: nodedof|=UX|UY|UZ; break;
        default: std::cerr << "Unknown element type" << std::endl; break;
      }
    }
    for (const auto node:element.second.nodes) {
      AnasysNodeDof[node] |= nodedof;
    }
  }



  std::vector<int> NodeIndex;
  N = 0;
  for (const auto &node:AnasysNodes) {
    int NodeDof=0;
    int numberOfLoaclDof =  __builtin_popcount(AnasysNodeDof[node.id]);
    NodeIndex.resize(numberOfLoaclDof);
    for (int i=0;i<numberOfLoaclDof;i++) {
      NodeDof = N + i + 1;
      NodeIndex[i] = NodeDof;
    }
    N+=numberOfLoaclDof;
    dof_manager.addNode(node.id, NodeIndex);
  }



  std::map<int, unsigned int> boundaryType;
  for (const auto &boundarytype : BoundaryAnasysData) {
    unsigned int boundaryDof = 0;
    if (boundarytype.dofLabel == "UX") {
      boundaryDof |=UX;
    } else if (boundarytype.dofLabel == "UY") {
      boundaryDof |=UY;
    } else if (boundarytype.dofLabel == "UZ") {
      boundaryDof |=UZ;
    } else if (boundarytype.dofLabel == "ROTX") {
      boundaryDof |=ROTX;
    } else if (boundarytype.dofLabel == "ROTY") {
      boundaryDof |=ROTY;
    } else if (boundarytype.dofLabel == "ROTZ") {
      boundaryDof |=ROTZ;
    } else {
      std::runtime_error("Nnkown boundary type");
    }
    boundaryType[boundarytype.nodeId] |= boundaryDof;
  }


  for (const auto &boundary:boundaryType) {
    if (boundary.second & UX) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][0]);
    }
   if (boundary.second & UY) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][1]);
    }
    if (boundary.second & UZ) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][2]);
    }
    if (boundary.second & ROTX) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][3]);
    }
    if (boundary.second & ROTY) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][4]);
    }
    if (boundary.second & ROTZ) {
      dof_manager.constrainDoF(dof_manager.node_to_dofs[boundary.first][5]);
    }
  }


  std::map<std::string, int> boundaryDirection;
  boundaryDirection["UX"] = 1;
  boundaryDirection["UY"] = 2;
  boundaryDirection["UZ"] = 3;
  boundaryDirection["ROTX"] = 4;
  boundaryDirection["ROTY"] = 5;
  boundaryDirection["ROTZ"] = 6;



  std::vector<int> ElementNodeIndex;

  for (const auto &element : elements) {

    ElementNodeIndex.reserve(AnasysElements[element.first]->nr);
    if (AnasysElementType.getElement(element.second.typ)->elementName=="COMBIN14") {
      for (const auto &node:element.second.nodes) {
        int keypoint=AnasysElementType.getElement(element.second.typ)->keyopts[1];
        if (keypoint==0) {
          ElementNodeIndex.push_back(dof_manager.node_to_dofs[node][0]);
          ElementNodeIndex.push_back(dof_manager.node_to_dofs[node][1]);
          ElementNodeIndex.push_back(dof_manager.node_to_dofs[node][2]);
        }
        else {
          ElementNodeIndex.push_back(dof_manager.node_to_dofs[node][keypoint-1]);
        }
      }
    }
    else {
      for (const auto &node : element.second.nodes) {
        std::vector<int> temp = dof_manager.node_to_dofs[node];
        for (size_t i = 0; i < AnasysElements[element.first]->dof; i++) {
          ElementNodeIndex.push_back(temp[i]);
        }
      }
    }

    dof_manager.addCell(element.first, ElementNodeIndex);
    ElementNodeIndex.clear();
  }

  LOG_INFO("Number of DoFs: {}" ,N);


  for (size_t i = 1; i <= N; i++) {
    dof_map.insert(i, i);
  }

  if (CE_enable == true) {
    for (const auto &CE : equations) {
      std::vector<int> slavenode;
      std::vector<double> CEvalue;
      int masterNode = 0;
      slavenode.clear();
      CEvalue.clear();
      for (const auto &term : CE.terms) {
        if (masterNode == 0) {
          masterNode = (old_to_new[term.node] - 1) *
                           dof_manager.node_to_dofs[term.node].size() +
                       boundaryDirection[term.direction];
        } else {
          slavenode.push_back((old_to_new[term.node] - 1) *
                                  dof_manager.node_to_dofs[term.node].size() +
                              boundaryDirection[term.direction]);
          CEvalue.push_back(-1.0 * term.coefficient);
        }
      }
      dof_manager.addCE(masterNode, slavenode, CEvalue);
    }
  }

  if (CP_enable) {
    for (const auto &set : coupledSets) {
      for (size_t index = 1; index < set.total_nodes; index++) {
        dof_map.updateByKey(
            (set.nodes[index] - 1) *
                    dof_manager.node_to_dofs[set.nodes[index]].size() +
                boundaryDirection[set.direction],
            (set.nodes[0] - 1) * dof_manager.node_to_dofs[set.nodes[0]].size() +
                boundaryDirection[set.direction]);
        dof_manager.cp_dofs_list.push_back(
            (set.nodes[index] - 1) *
                dof_manager.node_to_dofs[set.nodes[index]].size() +
            boundaryDirection[set.direction]);
        dof_manager.CP_association
            [(set.nodes[index] - 1) *
                 dof_manager.node_to_dofs[set.nodes[index]].size() +
             boundaryDirection[set.direction]] =
            (set.nodes[0] - 1) * dof_manager.node_to_dofs[set.nodes[0]].size() +
            boundaryDirection[set.direction];
      }
    }
  }
}


void AnsysLoad::saveToVTK() {
  vtkname=filename+".vtk";
  std::ofstream file(vtkname,std::ios::out);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  file<<"# vtk DataFile Version 3.0"<<std::endl;
  file<<"Generated by PMTOP"<<std::endl;
  file<<"ASCII"<<std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<"POINTS "<<AnasysNodes.size()<<" double"<<std::endl;
  for (const auto &node : AnasysNodes) {
    file<<node.x<<" "<<node.y<<" "<<node.z<<std::endl;
  }
  int numSize=0;
  for (const auto &element : AnasysElements) {
    numSize+=element.second->nodes.size()+1;
  }
  file<<"CELLS "<<AnasysElements.size()<<" "<<numSize<<std::endl;
  for (const auto &element : AnasysElements) {
    file<<element.second->nodes.size()<<" ";
    for (const auto &node : element.second->nodes) {
      file<<node.id-1<<" ";
    }
    file<<std::endl;
  }
  file<<"CELL_TYPES "<<AnasysElements.size()<<std::endl;
  for (const auto &element : AnasysElements) {
    file<<AnasysToVTKFormat(element.second->jtyp)<<std::endl;
  }
  file<<"CELL_DATA "<<AnasysElements.size()<<std::endl;
  file<<"SCALARS material int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.mat<<std::endl;
  }
  file<<"SCALARS ElementType int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.typ<<std::endl;
  }
  file<<"SCALARS Parts int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.sec<<std::endl;
  }
  file.close();
}

void AnsysLoad::saveToVTK(std::string name) {
  vtkname=filename+"_"+name+".vtk";
  std::ofstream file(vtkname,std::ios::out);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  file<<"# vtk DataFile Version 3.0"<<std::endl;
  file<<"Generated by PMTOP"<<std::endl;
  file<<"ASCII"<<std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<"POINTS "<<AnasysNodes.size()<<" double"<<std::endl;
  for (const auto &node : AnasysNodes) {
    file<<node.x<<" "<<node.y<<" "<<node.z<<std::endl;
  }
  int numSize=0;
  for (const auto &element : AnasysElements) {
    numSize+=element.second->nodes.size()+1;
  }
  file<<"CELLS "<<AnasysElements.size()<<" "<<numSize<<std::endl;
  for (const auto &element : AnasysElements) {
    file<<element.second->nodes.size()<<" ";
    for (const auto &node : element.second->nodes) {
      file<<node.id-1<<" ";
    }
    file<<std::endl;
  }
  file<<"CELL_TYPES "<<AnasysElements.size()<<std::endl;
  for (const auto &element : AnasysElements) {
    file<<AnasysToVTKFormat(element.second->jtyp)<<std::endl;
  }
  file<<"CELL_DATA "<<AnasysElements.size()<<std::endl;
  file<<"SCALARS material int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.mat<<std::endl;
  }
  file<<"SCALARS ElementType int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.typ<<std::endl;
  }
  file<<"SCALARS Parts int"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (const auto &element : elements) {
    file<<element.second.sec<<std::endl;
  }
  file.close();
}


void AnsysLoad::savePointValue(FloatArray& x,int Dim) {
  std::ofstream file(vtkname,std::ios::out | std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  file<<"POINT_DATA "<<AnasysNodes.size()<<std::endl;
  file<<"VECTORS displacement double"<<std::endl;
  for (int i=0;i<AnasysNodes.size();++i) {
    for (int j=0;j<Dim/2;j++) {
      file<<x[i*Dim+j]<<" ";
    }
    file<<std::endl;
  }
  file<<"VECTORS rotation double"<<std::endl;
  for (int i=0;i<AnasysNodes.size();++i) {
    for (int j=0;j<Dim/2;j++) {
      file<<x[i*Dim+j+Dim/2]<<" ";
    }
    file<<std::endl;
  }
  file.close();
}

void AnsysLoad::savePointValueShell(IntArray& index, FloatMatrix& vector) {
  std::ofstream file(vtkname,std::ios::out | std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  FloatArray SumE(NumberOfNode);
  SumE=0.0;
  file<<"POINT_DATA "<<AnasysNodes.size()<<std::endl;
  for (int ind=0;ind<index.giveSize();ind++) {
    file<<"VECTORS EigenVectorDisplacement"+std::to_string(index[ind])+" double"<<std::endl;
    for (int i=0;i<AnasysNodes.size();++i) {
      double SumL=0.0;
      for (int j=0;j<3;j++) {
        file<<vector(i*6+j,ind)<<" ";
        SumL+=vector(i*6+j,ind)*vector(i*6+j,ind);
      }
      SumE[i]+=std::sqrt(SumL);
      file<<std::endl;
    }
    file<<"VECTORS EigenVectorRotation"+std::to_string(index[ind])+" double"<<std::endl;
    for (int i=0;i<AnasysNodes.size();++i) {
      for (int j=0;j<3;j++) {
        file<<vector(i*6+j+3,ind)<<" ";
      }
      file<<std::endl;
    }
  }

  file<<"SCALARS EigenVectorSum double"<<std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for (int i=0;i<AnasysNodes.size();++i) {
    file<<SumE[i]<<" ";
    file<<std::endl;
  }

  file.close();
}

void AnsysLoad::savePointValueSolid(IntArray& index, FloatMatrix& vector) {
  std::ofstream file(vtkname,std::ios::out | std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  FloatArray SumE(NumberOfNode);
  SumE=0.0;
  file<<"POINT_DATA "<<AnasysNodes.size()<<std::endl;
  for (int ind=0;ind<index.giveSize();ind++) {
    file<<"VECTORS EigenVectorDisplacement"+std::to_string(index[ind])+" double"<<std::endl;
    for (int i=0;i<AnasysNodes.size();++i) {
      double SumL=0.0;
      for (int j=0;j<3;j++) {
        file<<vector(i*3+j,ind)<<" ";
        SumL+=vector(i*3+j,ind)*vector(i*3+j,ind);
      }
      SumE[i]+=std::sqrt(SumL);
      file<<std::endl;
    }
  }

  file<<"SCALARS EigenVectorSum double"<<std::endl;
  file << "LOOKUP_TABLE default" << std::endl;
  for (int i=0;i<AnasysNodes.size();++i) {
    file<<SumE[i]<<" ";
    file<<std::endl;
  }

  file.close();
}


void AnsysLoad::saveElementValue(std::string name,FloatArray& x) {
  std::ofstream file(vtkname,std::ios::out | std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  file<<"SCALARS "<<name<<" double"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (int i=0;i<AnasysElements.size();++i) {
    file<<x[i]<<std::endl;
  }
  file.close();
}

void AnsysLoad::saveElementValue(FloatArray& x) {
  std::ofstream file(vtkname,std::ios::out | std::ios::app);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << filename + ".vtk" << std::endl;
  }
  file<<"SCALARS ElementResult double"<<std::endl;
  file<<"LOOKUP_TABLE default"<<std::endl;
  for (int i=0;i<AnasysElements.size();++i) {
    file<<x[i]<<std::endl;
  }
  file.close();
}


IntArray AnsysLoad::getElementLocation(int i) {
  i=i+1;
  int size=AnasysElements[i]->nr;
  std::vector<Node> nodes = AnasysElements[i]->nodes;

  AnasysElements[i]->setMatrixLoacltion();


  std::vector<int> res_vector=AnasysElements[i]->MatrixLocaltion;
  IntArray res(size);
  for (size_t index = 0; index < size; index++) {
    res[index]=res_vector[index]-1;
  }
  return res;
}

void AnsysLoad::getConstrainedDofs(std::unordered_set<int> &ConstrainedDofs) {
  ConstrainedDofs.clear();
  for (const auto &[first, second] : dof_manager.constrained_dofs) {
    ConstrainedDofs.insert(first - 1);
  }
  if (getCP_enable()) {
    for (const auto &constraint : dof_manager.cp_dofs_list) {
      ConstrainedDofs.insert(constraint - 1);
    }
  }
}


void AnsysLoad::ReadAnasysElementStifAndElementMass() {
  std::ifstream infile(filename + "_elem_k.txt");
  if (!infile.is_open()) {
    std::cerr << "Error: Unable to open elem_k file!" << std::endl;
  }

  std::string line;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string word;
    int elem;
    int Index;
    if (line.find("Element") != std::string::npos) {

      iss >> word >> elem;
      Index = 0;
    } else {
      double value;
      iss >> value;

      AnasysElements[elem]->setStiff(Index, value);
      Index++;
    }
  }

  infile.close();

  std::ifstream infile2(filename + "_elem_m.txt");
  if (!infile2.is_open()) {
    std::cerr << "Error: Unable to open elem_m file!" << std::endl;
  }

  while (std::getline(infile2, line)) {
    std::istringstream iss(line);
    std::string word;
    int elem;
    int Index;
    if (line.find("Element") != std::string::npos) {
      iss >> word >> elem;
      Index = 0;
    } else {
      double value;
      iss >> value;
      AnasysElements[elem]->setMass(Index, value);
      Index++;
    }
  }

  infile2.close();
  LOG_INFO("Element Stif And Element Mass Reader Success!");

  std::ifstream infile3(filename + "_elem_c.txt");
  if (!infile3.is_open()) {
    std::cerr << "Error: Unable to open elem_c file!" << std::endl;
  }

  if (infile3.is_open()) {
    while (std::getline(infile3, line)) {
      std::istringstream iss(line);
      std::string word;
      int elem;
      int Index;
      if (line.find("Element") != std::string::npos) {
        iss >> word >> elem;
        Index = 0;
      } else {
        double value;
        iss >> value;
        AnasysElements[elem]->setDamp(Index, value);
        Index++;
      }
    }
    LOG_INFO("Element Damp Reader Success!");
  }

  infile3.close();

}

FloatMatrix AnsysLoad::getElementStif(int i, double x) {
  i=i+1;
  int size=AnasysElements[i]->nr;
  FloatMatrix res(size,size);
  for (size_t indexi = 0; indexi < size; indexi++) {
    for (size_t indexj = 0; indexj < size; indexj++) {
      res(indexi,indexj)=AnasysElements[i]->elementStif[indexj*size+indexi];
    }
  }
  res.times(x);
  return res;
}

FloatMatrix AnsysLoad::getElementMass(int i, double x) {
  i=i+1;
  int size=AnasysElements[i]->nr;
  FloatMatrix res(size,size);
  for (size_t indexi = 0; indexi < size; indexi++) {
    for (size_t indexj = 0; indexj < size; indexj++) {
      res(indexi,indexj)=AnasysElements[i]->elementMass[indexj*size+indexi];
    }
  }
  res.times(x);
  return res;
}

FloatMatrix AnsysLoad::getElementDamp(int i, double x) {
  i=i+1;
  int size=AnasysElements[i]->nr;
  FloatMatrix res(size,size);
  for (size_t indexi = 0; indexi < size; indexi++) {
    for (size_t indexj = 0; indexj < size; indexj++) {
      res(indexi,indexj)=AnasysElements[i]->elementDamp[indexj*size+indexi];
    }
  }
  res.times(x);
  return res;
}


KDNode* AnsysLoad::buildKDTree(std::vector<Node>& points, int dim, int depth) {
  if (points.empty()) return nullptr;

  std::vector<Node> pointsCopy = points;  // 拷贝一份，避免修改原数据
  const int k = dim;
  int axis = depth % k;

  int medianIndex = pointsCopy.size() / 2;
  std::nth_element(pointsCopy.begin(), pointsCopy.begin() + medianIndex, pointsCopy.end(),
                   Compare(axis));

  KDNode* node = new KDNode(pointsCopy[medianIndex]);

  std::vector<Node> leftPoints(pointsCopy.begin(), pointsCopy.begin() + medianIndex);
  std::vector<Node> rightPoints(pointsCopy.begin() + medianIndex + 1, pointsCopy.end());

  node->left = buildKDTree(leftPoints, dim, depth + 1);
  node->right = buildKDTree(rightPoints, dim, depth + 1);

  return node;
}



std::vector<Node> AnsysLoad::rangeSearch(KDNode *node, const Node &target, double radius, int depth, int dim) {
  std::vector<Node> result;

  // 如果当前节点为空，返回空的结果集
  if (node == nullptr)
    return result;

  const int k = dim;
  int axis = depth % k;   // 选择分割轴

  // 计算当前点到目标点的欧氏距离
  double dist = node->node.euclideanDistance(target);
  // 如果距离小于半径，则将当前节点加入结果
  if (dist <= radius) {
    result.push_back(node->node);
  }

  // 根据当前分割轴决定递归搜索左右子树
  if (axis == 0) { // x 轴
    if (target.x - radius < node->node.x) {
      auto leftResult = rangeSearch(node->left, target, radius, depth + 1, dim);
      result.insert(result.end(), leftResult.begin(), leftResult.end());
    }
    if (target.x + radius > node->node.x) {
      auto rightResult = rangeSearch(node->right, target, radius, depth + 1, dim);
      result.insert(result.end(), rightResult.begin(), rightResult.end());
    }
  }
  else if (axis == 1) { // y 轴
    if (target.y - radius < node->node.y) {
      auto leftResult = rangeSearch(node->left, target, radius, depth + 1, dim);
      result.insert(result.end(), leftResult.begin(), leftResult.end());
    }
    if (target.y + radius > node->node.y) {
      auto rightResult = rangeSearch(node->right, target, radius, depth + 1, dim);
      result.insert(result.end(), rightResult.begin(), rightResult.end());
    }
  }
  else if (axis == 2) { // z 轴
    if (target.z - radius < node->node.z) {
      auto leftResult = rangeSearch(node->left, target, radius, depth + 1, dim);
      result.insert(result.end(), leftResult.begin(), leftResult.end());
    }
    if (target.z + radius > node->node.z) {
      auto rightResult = rangeSearch(node->right, target, radius, depth + 1, dim);
      result.insert(result.end(), rightResult.begin(), rightResult.end());
    }
  }

  return result;
}

std::vector<Node> AnsysLoad::getCentroid() {
  std::vector<Node> res;
  for (const auto &elem : AnasysElements) {
    double idx=0.0;
    double dx=0.0,dy=0.0,dz=0.0,dthxy=0.0,dthyz=0.0,dthzx=0.0;
    for (const auto &node:elem.second->nodes) {
      idx+=1;
      dx+=node.x;
      dy+=node.y;
      dz+=node.z;
      dthxy+=node.thxy;
      dthyz+=node.thyz;
      dthzx+=node.thzx;
    }
    dx/=idx;dy/=idx;dz/=idx;dthxy/=idx;dthyz/=idx;dthzx/=idx;
    res.emplace_back(elem.first,dx,dy,dz,dthxy,dthyz,dthzx);
  }
  return res;
}

void AnsysLoad::printElementMatrix() {
  std::fstream f;
  f.open(filename + "_ElementMatrix.txt", std::ios::out);
  f << std::setprecision(20);
  for (const auto &element : AnasysElements) {
    if (element.second->MatrixLocaltion.empty()) {
      element.second->setMatrixLoacltion();
    }
    element.second->saveToStream(f);
  }

  f.close();

  std::cout << "Print ElementMatrix Success!" << std::endl;
}

bool ANSYSElementReader::parseElementLine(const std::string& line, ElementType& element) {
  // 查找 "ELEMENT TYPE" 关键字
  size_t typePos = line.find("ELEMENT TYPE");
  if (typePos == std::string::npos) return false;

  // 提取类型编号
  size_t numStart = typePos + 12; // "ELEMENT TYPE"后面的位置
  while (numStart < line.length() && !isdigit(line[numStart])) numStart++;

  std::string numStr = "";
  size_t pos = numStart;
  while (pos < line.length() && isdigit(line[pos])) {
    numStr += line[pos++];
  }

  if (!numStr.empty()) {
    element.typeNumber = std::stoi(numStr);
  }

  // 查找 "IS" 关键字后的元素名称和描述
  size_t isPos = line.find(" IS ");
  if (isPos != std::string::npos) {
    std::string remainder = line.substr(isPos + 4);
    std::istringstream iss(remainder);
    iss >> element.elementName;

    // 获取描述部分
    size_t nameEnd = remainder.find(element.elementName) + element.elementName.length();
    if (nameEnd < remainder.length()) {
      element.description = remainder.substr(nameEnd);
      // 去掉前后空白
      size_t start = element.description.find_first_not_of(" \t");
      size_t end = element.description.find_last_not_of(" \t\r\n");
      if (start != std::string::npos && end != std::string::npos) {
        element.description = element.description.substr(start, end - start + 1);
      }
    }
  }

  return true;
}

// 解析KEYOPT行
bool ANSYSElementReader::parseKeyoptLine(const std::string& line, int keyopts[], int startIndex) {
  // 查找等号后的数值
  size_t equalPos = line.find('=');
  if (equalPos == std::string::npos) return false;

  std::string values = line.substr(equalPos + 1);
  std::istringstream iss(values);

  int count = 0;
  int value;
  while (iss >> value && count < 6) {  // 每行最多6个值
    if (startIndex + count < 18) {
      keyopts[startIndex + count] = value;
    }
    count++;
  }

  return count > 0;
}

// 读取文件
bool ANSYSElementReader::readFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "无法打开文件: " << filename << std::endl;
    return false;
  }

  std::string line;
  ElementType currentElement;
  bool readingElement = false;
  int keyoptLineCount = 0;

  while (std::getline(file, line)) {
    // 跳过空行和注释行
    if (line.empty() || line.find("LIST ELEMENT TYPES") != std::string::npos ||
        line.find("CURRENT NODAL DOF") != std::string::npos ||
        line.find("THREE-DIMENSIONAL") != std::string::npos) {
      continue;
        }

    // 检查是否是元素类型行
    if (line.find("ELEMENT TYPE") != std::string::npos) {
      // 如果之前在读取元素，先保存
      if (readingElement) {
        elements.push_back(currentElement);
      }

      // 开始读取新元素
      currentElement = ElementType();
      if (parseElementLine(line, currentElement)) {
        readingElement = true;
        keyoptLineCount = 0;
      }
    }
    // 检查是否是KEYOPT行
    else if (line.find("KEYOPT") != std::string::npos && readingElement) {
      int startIndex = keyoptLineCount * 6;  // 每行6个值
      if (parseKeyoptLine(line, currentElement.keyopts, startIndex)) {
        keyoptLineCount++;
      }
    }
  }

  // 保存最后一个元素
  if (readingElement) {
    elements.push_back(currentElement);
  }

  file.close();
  return true;
}

// 打印所有元素信息
void ANSYSElementReader::printElements() const {
  std::cout << "读取到 " << elements.size() << " 个元素类型:\n" << std::endl;

  for (const auto& element : elements) {
    std::cout << "元素类型 " << element.typeNumber << ": "
              << element.elementName << " - " << element.description << std::endl;

    std::cout << "KEYOPT值:" << std::endl;
    for (int i = 0; i < 18; i++) {
      std::cout << std::setw(3) << element.keyopts[i];
      if ((i + 1) % 6 == 0) {  // 每6个值换行
        std::cout << std::endl;
      }
      else {
        std::cout << " ";
      }
    }
    std::cout << std::endl;
  }
}

// 获取元素数量
size_t ANSYSElementReader::getElementCount() const {
  return elements.size();
}

// 获取特定元素
const ElementType* ANSYSElementReader::getElement(int typeNumber) const {
  for (const auto& element : elements) {
    if (element.typeNumber == typeNumber) {
      return &element;
    }
  }
  return nullptr;
}

// 获取所有元素
const std::vector<ElementType>& ANSYSElementReader::getAllElements() const {
  return elements;
}