#include "PrismNedelec.h"

namespace JAUMIN {
namespace appu {
PrismNedelec::PrismNedelec(const hier::Patch<NDIM>& patch,
                           tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order,
                           tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord) :
  d_patch(patch),
  d_edge_order(edge_order),
  d_node_coord(node_coord)
{
  patch.getPatchTopology()->getCellAdjacencyEdges(d_cell_edge_ext, d_cell_edge_idx);
  patch.getPatchTopology()->getEdgeAdjacencyNodes(d_edge_node_ext, d_edge_node_idx);
  patch.getPatchTopology()->getCellAdjacencyNodes(d_cell_node_ext, d_cell_node_idx);
}
int PrismNedelec::nbas() const {
    return 9; // 三棱柱有 9 条边
}

int PrismNedelec::dim() const {
    return NDIM;
}
void PrismNedelec::basis(const int cell, const double* lambda, double* values) const
{
  double (*phi)[NDIM] = (double (*)[NDIM])values;

  double xi = lambda[0];
  double eta = lambda[1];
  double zeta = lambda[2]; // zeta in [-1, 1]

  // 1. 准备参考单元的重心坐标 L 和梯度
  double L[3] = {1.0 - xi - eta, xi, eta};
  double gradL[3][3] = {
    {-1.0, -1.0, 0.0},
    { 1.0,  0.0, 0.0},
    { 0.0,  1.0, 0.0}
  };

  // 2. 获取当前单元 6 个顶点的全局物理坐标
  double P[6][3];
  for(int n = 0; n < 6; n++) {
    int node_id = d_cell_node_idx[d_cell_node_ext[cell] + n];
    P[n][0] = (*d_node_coord)(0, node_id);
    P[n][1] = (*d_node_coord)(1, node_id);
    P[n][2] = (*d_node_coord)(2, node_id);
  }

  // 3. 计算 6 个节点的 H1 标量形函数的梯度 (用于求雅可比)
  double dN[6][3];
  for(int i = 0; i < 3; i++) {
    // 底面节点 (0,1,2): N = L_i * (1-zeta)/2
    dN[i][0] = gradL[i][0] * (1.0 - zeta) / 2.0;
    dN[i][1] = gradL[i][1] * (1.0 - zeta) / 2.0;
    dN[i][2] = -L[i] / 2.0;
    // 顶面节点 (3,4,5): N = L_i * (1+zeta)/2
    dN[i+3][0] = gradL[i][0] * (1.0 + zeta) / 2.0;
    dN[i+3][1] = gradL[i][1] * (1.0 + zeta) / 2.0;
    dN[i+3][2] = L[i] / 2.0;
  }

  // 4. 计算实时雅可比矩阵 J = \sum P_k \otimes \nabla N_k
  double J[3][3] = {0};
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 6; k++) {
        J[i][j] += P[k][i] * dN[k][j];
      }
    }
  }

  // 5. 计算 J 的行列式和逆矩阵 J^{-1}
  double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

  double invJ[3][3];
  invJ[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) / detJ;
  invJ[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2]) / detJ;
  invJ[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) / detJ;
  invJ[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2]) / detJ;
  invJ[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) / detJ;
  invJ[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2]) / detJ;
  invJ[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0]) / detJ;
  invJ[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1]) / detJ;
  invJ[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) / detJ;

  // 6. 为 9 条边计算基函数
  for(int ie = 0; ie < 9; ie++) {
    int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];

    // 方向处理
    int in_order = (*d_edge_order)(0, edge);
    if(!in_order) {
      int t = v0; v0 = v1; v1 = t;
    }

    // 寻找局部节点索引 [0, 5]
    int lv0 = -1, lv1 = -1;
    for(int n = 0; n < 6; n++) {
      int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
      if(nidx == v0) lv0 = n;
      if(nidx == v1) lv1 = n;
    }

    double E_ref[3] = {0, 0, 0};

    // 判断属于哪种类型的边，计算参考空间的 E_ref
    if(lv0 < 3 && lv1 < 3) {
      // 类型 1: 底面水平边
      double z_fac = (1.0 - zeta) / 2.0;
      for(int d=0; d<3; d++) {
        E_ref[d] = (L[lv0] * gradL[lv1][d] - L[lv1] * gradL[lv0][d]) * z_fac;
      }
    }
    else if(lv0 >= 3 && lv1 >= 3) {
      // 类型 2: 顶面水平边
      int a = lv0 - 3;
      int b = lv1 - 3;
      double z_fac = (1.0 + zeta) / 2.0;
      for(int d=0; d<3; d++) {
        E_ref[d] = (L[a] * gradL[b][d] - L[b] * gradL[a][d]) * z_fac;
      }
    }
    else {
      // 类型 3: 垂直侧边 (nabla zeta = [0,0,1]^T)
      if(lv1 == lv0 + 3) { // 向上
        E_ref[2] = L[lv0] * 0.5;
      } else if(lv0 == lv1 + 3) { // 向下
        E_ref[2] = -L[lv1] * 0.5;
      }
    }

    // 7. 使用逆雅可比矩阵的转置进行协变 Piola 变换: E_phys = J^{-T} * E_ref
    phi[ie][0] = invJ[0][0]*E_ref[0] + invJ[1][0]*E_ref[1] + invJ[2][0]*E_ref[2];
    phi[ie][1] = invJ[0][1]*E_ref[0] + invJ[1][1]*E_ref[1] + invJ[2][1]*E_ref[2];
    phi[ie][2] = invJ[0][2]*E_ref[0] + invJ[1][2]*E_ref[1] + invJ[2][2]*E_ref[2];
  }
} // end of basis function
void PrismNedelec::curl(const int cell, const double* lambda, double* values) const{

}
}
}
