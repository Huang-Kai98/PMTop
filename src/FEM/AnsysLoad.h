#pragma once

#include <fstream>
#include <intarray.h>
#include "floatmatrix.h"
#include "Boundary.h"
#include "Node.h"
#include "DofManager.h"
#include "Element.h"
#include <memory>
#include "bimap.h"
#include <map>
#include <set>

#include "Config.h"
#include "unordered_set"
#include "Logger.h"

struct ElementNodesInfo {
  int elem;
  int mat;
  int typ;
  int rel;
  int esy;
  int sec;
  std::vector<int> nodes;
};

enum NodeDof {
  UX      = 1 << 0,
  UY = 1 << 1,
  UZ    = 1 << 2,
  ROTX     = 1 << 3,
  ROTY     = 1 << 4,
  ROTZ     = 1 << 5,
};

struct Term {
  int node;
  std::string direction;
  double coefficient;
};

struct ConstraintEquation {
  int equation_number;
  int term_count;
  double constant;
  std::vector<Term> terms;
};
struct CoupledSet {
  int set_id;
  std::string direction;
  int total_nodes;
  std::vector<int> nodes;
};

struct KDNode {
  Node node;             // d 维坐标
  KDNode* left;          // 左子树
  KDNode* right;         // 右子树

  // 构造函数
  KDNode(const Node& pt) : node(pt), left(nullptr), right(nullptr) {}
};

// 比较函数：用于在指定维度上对点进行排序
struct Compare {
  int axis;
  Compare(int axis) : axis(axis) {}

  bool operator()(const Node& a, const Node& b) {
    if (axis == 0) return a.x < b.x;
    if (axis == 1) return a.y < b.y;
    return a.z < b.z;
  }
};

struct ElementType {
    int typeNumber;
    std::string elementName;
    std::string description;
    int keyopts[18];  // 18个KEYOPT值

    ElementType() {
        typeNumber = 0;
        // 初始化keyopts数组
        for(int i = 0; i < 18; i++) {
            keyopts[i] = 0;
        }
    }
};
// 1. 结构体 (保持精简，只存必要数据)
struct ElementData {
    int nr = 0;             // 从文件中读取并存储
    std::vector<double> K;  // 刚度
    std::vector<double> M;  // 质量
    std::vector<double> C;  // 阻尼
};

// 2. 数据库
class MatrixDatabase {
    inline static std::map<int, ElementData> store;
public:
    static ElementData& getEntry(int tag) {
        return store[tag];
    }
    static bool exists(int tag) { return store.count(tag); }
};

// 3. 加载器 (重新加入 NR 读取逻辑)
class MatrixLoader {
public:
    static bool load(const std::string& filepath, int targetTag = -1) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "File not found: " << filepath << std::endl;
            return false;
        }

        std::string token;
        int currentTag = -1;
        bool skip = false;
        int currentMatrixSize = 0; // 动态计算

        while (file >> token) {
            // --- 1. 处理注释 ---
            if (token[0] == '#') {
                std::string dummy; std::getline(file, dummy); continue;
            }

            // --- 2. 处理 TAG ---
            if (token == "TAG") {
                file >> currentTag;
                currentMatrixSize = 0; // 切换Tag时重置大小，防止沿用上一个Tag的数据

                // 决定是否跳过 (如果是目标Tag或者没指定目标，就不跳过)
                skip = (targetTag != -1 && currentTag != targetTag);
            }
            else if (token == "NR") {
                int val;
                file >> val;

                // 即使 skip==true，也要计算大小，以便后续正确跳过矩阵数据
                currentMatrixSize = val * val;

                // 如果不跳过，顺便把 nr 存进数据库
                if (!skip) {
                    MatrixDatabase::getEntry(currentTag).nr = val;
                }
            }
            else if (token == "MATRIX") {
                std::string type;
                file >> type; // STIFFNESS, MASS, DAMP

                // 安全检查：如果还没读到 NR 就遇到 MATRIX，或者 NR=0
                if (currentMatrixSize == 0) {
                    std::cerr << "Error: MATRIX found before NR defined for TAG " << currentTag << std::endl;
                    return false;
                }

                // 如果处于跳过模式，仅读取并丢弃数据
                if (skip) {
                    double d;
                    for (int i = 0; i < currentMatrixSize; ++i) file >> d;
                    continue;
                }

                // 正常读取模式
                ElementData& data = MatrixDatabase::getEntry(currentTag);
                std::vector<double>* vec = nullptr;

                if (type == "STIFFNESS") vec = &data.K;
                else if (type == "MASS")     vec = &data.M;
                else if (type == "DAMP")     vec = &data.C;
                else {
                    // 未知类型，跳过
                    double d;
                    for(int i=0; i<currentMatrixSize; ++i) file >> d;
                    continue;
                }

                // 读取数据
                vec->resize(currentMatrixSize);
                for (int i = 0; i < currentMatrixSize; ++i) {
                    file >> (*vec)[i];
                }
            }
        }
        return true;
    }
};

class ANSYSElementReader {
private:
    std::vector<ElementType> elements;

    // 解析元素类型行
    bool parseElementLine(const std::string& line, ElementType& element);

    // 解析KEYOPT行
    bool parseKeyoptLine(const std::string& line, int keyopts[], int startIndex);

public:
    // 读取文件
    bool readFile(const std::string& filename) ;

    // 打印所有元素信息
    void printElements() const ;

    // 获取元素数量
    size_t getElementCount() const ;

    // 获取特定元素
    const ElementType* getElement(int typeNumber) const ;

    // 获取所有元素
    const std::vector<ElementType>& getAllElements() const ;
};

struct RayleighDampingCoefficients {
    double alpha;
    double beta;
};

struct ModelSquareSize {
    double length;
    double width;
};

class AnsysLoad {
private:
  std::string filename;
  //std::map<int, std::string> AnasysElementType;     // 存储单元类型信息
  ANSYSElementReader AnasysElementType;  // 存储单元类型信息
  std::map<int, std::unique_ptr<Element>> AnasysElements;
  std::vector<Node> AnasysNodes;                    // 存储节点信息
  std::map<int, unsigned int> AnasysNodeDof;
  int NumberOfNode;                                 // 节点数量
  std::map<int, ElementNodesInfo> elements;         // 存储单元节点信息
  int NumberOfElement;                              // 单元数量
  std::vector<std::string> BoundaryAnasysdofLabels; // 存储 DOF 标签
  std::vector<BoundaryData> BoundaryAnasysData;     // 存储边界条件信息
  std::vector<ConstraintEquation> equations;        // 存储 CE 单元信息
  std::vector<CoupledSet> coupledSets;

  bool CE_enable;
  bool CP_enable;
  size_t vertex_count; // 顶点数
  size_t edge_count;   // 边数

  std::map<int, int> old_to_new;
  std::map<int, int> new_to_old;
  std::map<int, int> MatrixMap;
  std::map<int, int> nodeDof;
  BiMap<int, int> dof_map;
  std::string vtkname;
  int N;
  friend class Top;
  friend class Sensitivity;
  friend class TopProblem;
public:
  DofManager dof_manager;
  std::map<int, std::set<int>> adjacency_list; // 节点邻接表
  std::map<int, std::set<int>> adjacent_cells; // 单元邻接表
  AnsysLoad(const Config& cfg) {
    // 1. 获取读取模式，默认为全量读取
    std::string read_mode = cfg.get<std::string>("read_mode");

    if (read_mode == "all_element") {
      this->filename = cfg.get<std::string>("filename");
      initFromFile(); // 调用全量初始化
    }
    else if (read_mode == "single_element") {
      int element_tag = cfg.get<int>("element_tag");
      this->filename = cfg.get<std::string>("filename"); // 标记一下，防止文件读取函数误用
      initSingleElement(element_tag); // 调用单单元初始化
    }
    else if (read_mode == "single_element_two_coeff") {
        int element_tag = cfg.get<int>("element_tag");
        this->filename = cfg.get<std::string>("filename"); // 标记一下，防止文件读取函数误用
        RayleighDampingCoefficients Coef1,Coef2;
        Coef1.alpha = cfg.get<double>("RayleighDampingCoef.Coef1.alpha");
        Coef1.beta = cfg.get<double>("RayleighDampingCoef.Coef1.beta");

        Coef2.alpha =cfg.get<double>("RayleighDampingCoef.Coef2.alpha");
        Coef2.beta = cfg.get<double>("RayleighDampingCoef.Coef2.beta");

        if (cfg.get<std::string>("Interface.Type")=="OneWidth") {
            auto Interface=cfg.get<double>("Interface.Interface");
            initSingleElement(element_tag,Coef1,Coef2,Interface);
        }
        else {
            LOG_ERROR("Nuknow Type!");
        }


    }
    else {
      LOG_ERROR("Unknown read_mode: {}",read_mode);
      // 可以选择抛出异常或者设置默认行为
    }
  };

    void initSingleElement(int element_tag,RayleighDampingCoefficients Coef1,RayleighDampingCoefficients Coef2,
        double Interface) {
        ReadAnsysNode();
        ReadAnsysElementType();
        ReadAnasysElementNodes();
        ReadAnasysBoundaryCondition();
        ReadCEelement();
        ReadCoupledSets();
        Adjacency_list_element();
        Adjacency_list_node();
        for (size_t i = 1; i <= getNumberOfNodes(); i++) {
            old_to_new[i] = i;
            new_to_old[i] = i;
        }
        dofManage();
        auto Centroid=getCentroid();

        LOG_INFO("Initializing single element mode for tag: {}" , std::to_string(element_tag));


        MatrixLoader::load("element_library.txt",element_tag);
        auto test=MatrixDatabase::exists(element_tag);
        auto ElementData=MatrixDatabase::getEntry(element_tag);
        auto ElementStif=ElementData.K;
        auto ElementMass=ElementData.M;
#pragma omp parallel for
        for (int i = 1; i <= getNumberOfElements(); i++) {



            for (int j = 0; j < AnasysElements[i]->nr; j++) {
                for (int k=0;k<AnasysElements[i]->nr;k++) {
                    AnasysElements[i]->setStiff(j,k,ElementStif[k*AnasysElements[i]->nr+j]);
                }
            }

            for (int j = 0; j < AnasysElements[i]->nr; j++) {
                for (int k=0;k<AnasysElements[i]->nr;k++) {
                    AnasysElements[i]->setMass(j,k,ElementMass[k*AnasysElements[i]->nr+j]);
                }
            }

            if (Centroid[i-1].y<=Interface) {
                elements[i].sec=1;
                for (int j = 0; j < AnasysElements[i]->nr; j++) {
                    for (int k=0;k<AnasysElements[i]->nr;k++) {

                        AnasysElements[i]->setDamp(j,k,ElementStif[k*AnasysElements[i]->nr+j]*Coef1.beta+
                            ElementMass[k*AnasysElements[i]->nr+j]*Coef1.alpha);
                    }
                }
            }
            else {
                elements[i].sec=2;
                for (int j = 0; j < AnasysElements[i]->nr; j++) {
                    for (int k=0;k<AnasysElements[i]->nr;k++) {

                        AnasysElements[i]->setDamp(j,k,ElementStif[k*AnasysElements[i]->nr+j]*Coef2.beta+
                            ElementMass[k*AnasysElements[i]->nr+j]*Coef2.alpha);
                    }
                }
            }
        }
        LOG_INFO("single model loaded successfully." );

    }


  // 处理单个单元的初始化逻辑
  void initSingleElement(int element_tag) {
    ReadAnsysNode();
    ReadAnsysElementType();
    ReadAnasysElementNodes();
    ReadAnasysBoundaryCondition();
    ReadCEelement();
    ReadCoupledSets();
    Adjacency_list_element();
    Adjacency_list_node();
    for (size_t i = 1; i <= getNumberOfNodes(); i++) {
      old_to_new[i] = i;
      new_to_old[i] = i;
    }
    dofManage();

    LOG_INFO("Initializing single element mode for tag: {}" , std::to_string(element_tag));


      MatrixLoader::load("element_library.txt",element_tag);
      auto test=MatrixDatabase::exists(element_tag);
      auto ElementData=MatrixDatabase::getEntry(element_tag);
        auto ElementStif=ElementData.K;
        auto ElementMass=ElementData.M;
        auto ElementDamp=ElementData.C;
#pragma omp parallel for
      for (int i = 1; i <= getNumberOfElements(); i++) {

          for (int j = 0; j < AnasysElements[i]->nr; j++) {
              for (int k=0;k<AnasysElements[i]->nr;k++) {
                  AnasysElements[i]->setStiff(j,k,ElementStif[k*AnasysElements[i]->nr+j]);
              }
          }

          for (int j = 0; j < AnasysElements[i]->nr; j++) {
              for (int k=0;k<AnasysElements[i]->nr;k++) {
                  AnasysElements[i]->setMass(j,k,ElementMass[k*AnasysElements[i]->nr+j]);
              }
          }

          for (int j = 0; j < AnasysElements[i]->nr; j++) {
              for (int k=0;k<AnasysElements[i]->nr;k++) {
                  AnasysElements[i]->setDamp(j,k,ElementDamp[k*AnasysElements[i]->nr+j]);
              }
          }
      }
      LOG_INFO("single model loaded successfully." );


  }
  void initFromFile() {
    ReadAnsysNode();
    ReadAnsysElementType();
    ReadAnasysElementNodes();
    ReadAnasysBoundaryCondition();
    ReadCEelement();
    ReadCoupledSets();
    Adjacency_list_element();
    Adjacency_list_node();
    for (size_t i = 1; i <= getNumberOfNodes(); i++) {
      old_to_new[i] = i;
      new_to_old[i] = i;
    }
    dofManage();
    ReadAnasysElementStifAndElementMass();
    LOG_INFO("Full model loaded successfully.");
  }
    std::string getElementType(int i){return AnasysElements[i+1]->jtyp;};

  ~AnsysLoad(){};
  ElementNodesInfo parseLine(const std::string &line);
  void ReadAnasysElementStifAndElementMass();
  void ReadAnsysElementType();        // 读取单元类型ETLIST
  void ReadAnsysNode();               // 读取节点信息NLIST
  void ReadAnasysElementNodes();      // 读取单元节点信息ELIST
  void ReadAnasysBoundaryCondition(); // 读取边界条件信息
  void ReadCEelement();               // 读取 CE 单元信息
  void ReadCoupledSets();              // 读取 CP 单元信息
  void Adjacency_list_node();         // 生成节点邻接表
  void Adjacency_list_element();      // 生成单元邻接表(共边)
  void Adjacency_list_element_node(); // 生成节点邻接表(共点)
  void dofManage();
  void saveToVTK();
  void saveToVTK(std::string name);
  void savePointValue(FloatArray& x,int Dim=6);
  void savePointValueShell(IntArray &index, FloatMatrix &vector);
  void savePointValueSolid(IntArray &index, FloatMatrix &vector);
  void saveElementValue(std::string name,FloatArray& x);
  void saveElementValue(FloatArray& x);
  static KDNode* buildKDTree(std::vector<Node>& points, int dim, int depth = 0);
  std::vector<Node> rangeSearch(KDNode *node,const Node &target, double radius, int depth, int dim = 3);
  // 递归打印 KD 树（前序遍历）
  void printKDTree(const KDNode* node,int level = 0) {
    if (node == nullptr) return;
    std::cout << std::string(level * 2, ' ') << "Level " << level << ": (";
    std::cout << node->node.id;
    std::cout << ")" << std::endl;
    printKDTree(node->left, level + 1);
    printKDTree(node->right, level + 1);
  }

  std::vector<Node> &getAnasysNodes() { return AnasysNodes; }
  std::map<int, ElementNodesInfo> getElements() { return elements; }
  std::vector<BoundaryData> getBoundaryAnasysData() {
    return BoundaryAnasysData;
  }
  std::vector<Node> getCentroid();
  std::vector<ConstraintEquation> getEquations() { return equations; }
  std::vector<CoupledSet> getCoupledSets() { return coupledSets; }

  std::map<int, std::set<int>> getAdjacencyList() { return adjacency_list; }
  std::map<int, std::set<int>> getAdjacentCells() { return adjacent_cells; }
  void getConstrainedDofs(std::unordered_set<int>& ConstrainedDofs);
  size_t getVertexCount() const { return vertex_count; }
  size_t getEdgeCount() const { return edge_count; }
  size_t getNumberOfNodes() const { return NumberOfNode; }
  size_t getNumberOfElements() const { return NumberOfElement; }
  size_t getNumberOfDof() const {return N;}
  bool getCE_enable() { return CE_enable; }
  bool getCP_enable() { return CP_enable; }
  IntArray getElementLocation(int i);
  FloatMatrix getElementStif(int i, double x);
  FloatMatrix getElementMass(int i, double x);
  FloatMatrix getElementDamp(int i, double x);
  static std::unique_ptr<Element> createElement(const std::string &type);
  void printElementMatrix();
};
