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

#ifdef DEBUG_CHECK_ASSERTIONS
  // 对第一个单元抽样校验节点排布约定:
  //   底面 = 局部索引 0,1,2 (逆时针)
  //   顶面 = 局部索引 3,4,5 (V3 在 V0 正上方, 以此类推)
  //   垂直边必须连接 V_i ↔ V_{i+3}
  if (d_cell_edge_ext.size() > 1) {
    const int cell = 0;
    for (int ie = 0; ie < NBAS; ie++) {
      const int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
      const int n0   = d_edge_node_idx[d_edge_node_ext[edge] + 0];
      const int n1   = d_edge_node_idx[d_edge_node_ext[edge] + 1];

      int l0 = -1, l1 = -1;
      for (int n = 0; n < 6; n++) {
        const int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
        if (nidx == n0) l0 = n;
        if (nidx == n1) l1 = n;
      }
      TBOX_ASSERT(l0 != -1 && l1 != -1);

      if (!(l0 < 3 && l1 < 3) && !(l0 >= 3 && l1 >= 3)) {
        // 垂直边: 必须满足 l1 == l0 + 3 或 l0 == l1 + 3
        TBOX_ASSERT(l1 == l0 + 3 || l0 == l1 + 3);
      }
    }
  }
#endif
}

int PrismNedelec::nbas() const {
    return NBAS;
}

int PrismNedelec::dim() const {
    return NDIM;
}

// ==========================================================================
// computeJacobian: basis() 与 curl() 共用的几何映射计算
// 在指定积分点计算参考空间到物理空间的 Jacobian、逆矩阵及行列式。
// ==========================================================================
void PrismNedelec::computeJacobian(
    const int cell, const double* lambda,
    double P[6][3], PrismJacobian& jac) const
{
  const double xi   = lambda[0];
  const double eta  = lambda[1];
  const double zeta = lambda[2];

  const double L[3] = {1.0 - xi - eta, xi, eta};

  // 三角形重心坐标的常数梯度 (参考空间)
  static const double gradL[3][3] = {
    {-1.0, -1.0, 0.0},
    { 1.0,  0.0, 0.0},
    { 0.0,  1.0, 0.0}
  };

  // 1. 获取当前单元 6 个顶点的全局物理坐标
  for (int n = 0; n < 6; n++) {
    int node_id = d_cell_node_idx[d_cell_node_ext[cell] + n];
    P[n][0] = (*d_node_coord)(0, node_id);
    P[n][1] = (*d_node_coord)(1, node_id);
    P[n][2] = (*d_node_coord)(2, node_id);
  }

  // 2. 计算 6 个节点 H1 标量形函数的梯度 dN_k / dξ
  //    底面 N_i   = L_i * (1-zeta)/2
  //    顶面 N_i+3 = L_i * (1+zeta)/2
  double dN[6][3];
  for (int i = 0; i < 3; i++) {
    dN[i][0]   = gradL[i][0] * (1.0 - zeta) * 0.5;
    dN[i][1]   = gradL[i][1] * (1.0 - zeta) * 0.5;
    dN[i][2]   = -L[i] * 0.5;
    dN[i+3][0] = gradL[i][0] * (1.0 + zeta) * 0.5;
    dN[i+3][1] = gradL[i][1] * (1.0 + zeta) * 0.5;
    dN[i+3][2] =  L[i] * 0.5;
  }

  // 3. 组装 Jacobian: J_{mn} = Σ_k P_k^m · dN_k^n
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      jac.J[i][j] = 0.0;
      for (int k = 0; k < 6; k++) {
        jac.J[i][j] += P[k][i] * dN[k][j];
      }
    }
  }

  // 4. 行列式
  jac.detJ = jac.J[0][0] * (jac.J[1][1]*jac.J[2][2] - jac.J[1][2]*jac.J[2][1])
           - jac.J[0][1] * (jac.J[1][0]*jac.J[2][2] - jac.J[1][2]*jac.J[2][0])
           + jac.J[0][2] * (jac.J[1][0]*jac.J[2][1] - jac.J[1][1]*jac.J[2][0]);

  const double invDet = 1.0 / jac.detJ;

  // 5. 逆矩阵 (手动 Cramer 法则, GCC 4.8.5 兼容)
  jac.invJ[0][0] = (jac.J[1][1]*jac.J[2][2] - jac.J[1][2]*jac.J[2][1]) * invDet;
  jac.invJ[0][1] = (jac.J[0][2]*jac.J[2][1] - jac.J[0][1]*jac.J[2][2]) * invDet;
  jac.invJ[0][2] = (jac.J[0][1]*jac.J[1][2] - jac.J[0][2]*jac.J[1][1]) * invDet;
  jac.invJ[1][0] = (jac.J[1][2]*jac.J[2][0] - jac.J[1][0]*jac.J[2][2]) * invDet;
  jac.invJ[1][1] = (jac.J[0][0]*jac.J[2][2] - jac.J[0][2]*jac.J[2][0]) * invDet;
  jac.invJ[1][2] = (jac.J[0][2]*jac.J[1][0] - jac.J[0][0]*jac.J[1][2]) * invDet;
  jac.invJ[2][0] = (jac.J[1][0]*jac.J[2][1] - jac.J[1][1]*jac.J[2][0]) * invDet;
  jac.invJ[2][1] = (jac.J[0][1]*jac.J[2][0] - jac.J[0][0]*jac.J[2][1]) * invDet;
  jac.invJ[2][2] = (jac.J[0][0]*jac.J[1][1] - jac.J[0][1]*jac.J[1][0]) * invDet;
}

// ==========================================================================
// basis: 计算 9 条边在指定积分点的 Nedelec 基函数值
//   E_ref → Piola 协变映射 (J^{-T} · E_ref) → E_phys
// ==========================================================================
void PrismNedelec::basis(const int cell, const double* lambda, double* values) const
{
  double (*phi)[NDIM] = (double (*)[NDIM])values;

  const double xi   = lambda[0];
  const double eta  = lambda[1];
  const double zeta = lambda[2];

  const double L[3] = {1.0 - xi - eta, xi, eta};

  static const double gradL[3][3] = {
    {-1.0, -1.0, 0.0},
    { 1.0,  0.0, 0.0},
    { 0.0,  1.0, 0.0}
  };

  // ★ 调用共用的 Jacobian 计算
  PrismJacobian jac;
  double P[6][3];
  computeJacobian(cell, lambda, P, jac);

  for (int ie = 0; ie < NBAS; ie++) {
    int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];

    // 边方向处理
    int in_order = (*d_edge_order)(0, edge);
    if (!in_order) {
      int t = v0; v0 = v1; v1 = t;
    }

    // 寻找局部节点索引 [0, 5]
    int lv0 = -1, lv1 = -1;
    for (int n = 0; n < 6; n++) {
      int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
      if (nidx == v0) lv0 = n;
      if (nidx == v1) lv1 = n;
    }

    double E_ref[3] = {0, 0, 0};

    if (lv0 < 3 && lv1 < 3) {
      // 类型 1: 底面水平边
      double z_fac = (1.0 - zeta) * 0.5;
      for (int d = 0; d < 3; d++) {
        E_ref[d] = (L[lv0] * gradL[lv1][d] - L[lv1] * gradL[lv0][d]) * z_fac;
      }
    }
    else if (lv0 >= 3 && lv1 >= 3) {
      // 类型 2: 顶面水平边
      int a = lv0 - 3;
      int b = lv1 - 3;
      double z_fac = (1.0 + zeta) * 0.5;
      for (int d = 0; d < 3; d++) {
        E_ref[d] = (L[a] * gradL[b][d] - L[b] * gradL[a][d]) * z_fac;
      }
    }
    else {
      // 类型 3: 垂直侧边
      if (lv1 == lv0 + 3) {
        E_ref[2] = L[lv0] * 0.5;
      } else if (lv0 == lv1 + 3) {
        E_ref[2] = -L[lv1] * 0.5;
      }
    }

    // Piola 协变映射: E_phys = J^{-T} · E_ref
    phi[ie][0] = jac.invJ[0][0]*E_ref[0] + jac.invJ[1][0]*E_ref[1] + jac.invJ[2][0]*E_ref[2];
    phi[ie][1] = jac.invJ[0][1]*E_ref[0] + jac.invJ[1][1]*E_ref[1] + jac.invJ[2][1]*E_ref[2];
    phi[ie][2] = jac.invJ[0][2]*E_ref[0] + jac.invJ[1][2]*E_ref[1] + jac.invJ[2][2]*E_ref[2];
  }
}

// ==========================================================================
// curl: 计算 9 条边在指定积分点的 Nedelec 基函数旋度
//   curl_ξ(E_ref) → Piola 逆变映射 ((1/detJ) · J · curl_ref) → curl_x(E_phys)
// ==========================================================================
void PrismNedelec::curl(const int cell, const double* lambda, double* values) const
{
  double (*curlphi)[NDIM] = (double (*)[NDIM])values;

  const double xi   = lambda[0];
  const double eta  = lambda[1];
  const double zeta = lambda[2];

  const double L[3] = {1.0 - xi - eta, xi, eta};

  static const double gradL[3][3] = {
    {-1.0, -1.0, 0.0},
    { 1.0,  0.0, 0.0},
    { 0.0,  1.0, 0.0}
  };

  // ★ 调用共用的 Jacobian 计算 (与 basis 完全一致)
  PrismJacobian jac;
  double P[6][3];
  computeJacobian(cell, lambda, P, jac);

  const double invDet = 1.0 / jac.detJ;

  for (int ie = 0; ie < NBAS; ie++) {
    int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];

    int in_order = (*d_edge_order)(0, edge);
    if (!in_order) {
      int t = v0; v0 = v1; v1 = t;
    }

    int lv0 = -1, lv1 = -1;
    for (int n = 0; n < 6; n++) {
      int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
      if (nidx == v0) lv0 = n;
      if (nidx == v1) lv1 = n;
    }

    double curl_E_ref[3] = {0, 0, 0};

    if (lv0 < 3 && lv1 < 3) {
      // 底面水平边: curl(B·A) = ∇B × A + B · curl(A)
      //   B = (1-zeta)/2, A = L_a·∇L_b - L_b·∇L_a, curl(A) = 2·∇L_a×∇L_b
      double A[3];
      for (int d = 0; d < 3; d++)
        A[d] = L[lv0] * gradL[lv1][d] - L[lv1] * gradL[lv0][d];

      double curlA[3];
      CROSS_PRODUCT(gradL[lv0], gradL[lv1], curlA);
      curlA[0] *= 2.0; curlA[1] *= 2.0; curlA[2] *= 2.0;

      double gradB[3] = {0, 0, -0.5};
      double term1[3];
      CROSS_PRODUCT(gradB, A, term1);

      double z_fac = (1.0 - zeta) * 0.5;
      for (int d = 0; d < 3; d++)
        curl_E_ref[d] = term1[d] + z_fac * curlA[d];
    }
    else if (lv0 >= 3 && lv1 >= 3) {
      // 顶面水平边: ∇B = (0, 0, +0.5)
      int a = lv0 - 3, b = lv1 - 3;
      double A[3];
      for (int d = 0; d < 3; d++)
        A[d] = L[a] * gradL[b][d] - L[b] * gradL[a][d];

      double curlA[3];
      CROSS_PRODUCT(gradL[a], gradL[b], curlA);
      curlA[0] *= 2.0; curlA[1] *= 2.0; curlA[2] *= 2.0;

      double gradB[3] = {0, 0, 0.5};
      double term1[3];
      CROSS_PRODUCT(gradB, A, term1);

      double z_fac = (1.0 + zeta) * 0.5;
      for (int d = 0; d < 3; d++)
        curl_E_ref[d] = term1[d] + z_fac * curlA[d];
    }
    else {
      // 垂直侧边: curl(L_i·0.5·e_z) = ∇L_i × (0.5·e_z)
      if (lv1 == lv0 + 3) {
        double V[3] = {0, 0, 0.5};
        CROSS_PRODUCT(gradL[lv0], V, curl_E_ref);
      } else if (lv0 == lv1 + 3) {
        double V[3] = {0, 0, -0.5};
        CROSS_PRODUCT(gradL[lv1], V, curl_E_ref);
      }
    }

    // Piola 逆变映射: curl_phys = (1/detJ) · J · curl_ref
    curlphi[ie][0] = invDet * (jac.J[0][0]*curl_E_ref[0] + jac.J[0][1]*curl_E_ref[1] + jac.J[0][2]*curl_E_ref[2]);
    curlphi[ie][1] = invDet * (jac.J[1][0]*curl_E_ref[0] + jac.J[1][1]*curl_E_ref[1] + jac.J[1][2]*curl_E_ref[2]);
    curlphi[ie][2] = invDet * (jac.J[2][0]*curl_E_ref[0] + jac.J[2][1]*curl_E_ref[1] + jac.J[2][2]*curl_E_ref[2]);
  }
}

} // namespace appu
} // namespace JAUMIN
