//
// 文件名:      PrismNedelec.h
// 描述  :      三棱柱 (Prism) 零阶 Nedelec 边缘元基函数
//

#ifndef included_appu_PrismNedelec
#define included_appu_PrismNedelec

#if NDIM == 2
#error FOR 3D CASE ONLY
#endif

#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "NodeData.h"
#include "JAUMIN_Macros.h"

namespace JAUMIN {
namespace appu {
class PrismNedelec
{
public:
    enum { NBAS = 9 };   // ★ GCC 4.8.5 安全的编译期常量，对齐四面体 Nedelec 的范式

    /**
         * 构造函数
         * 【改动】移除了 cell_jacobian，加入了 node_coord 用于实时计算双线性 Jacobian
         */
    PrismNedelec(const hier::Patch<NDIM>& patch,
                 tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order,
                 tbox::Pointer<pdat::NodeData<NDIM, double> > node_coord);

    /**
         * 单元内基函数的数量 (三棱柱有 9 条边)
         */
    int nbas() const;

    /**
         * 基函数的维数 (3D 矢量)
         */
    int dim() const;
    void basis(const int cell, const double* lambda, double* values) const;

    /**
         * 计算三棱柱 Nedelec 基函数的旋度 (需经过 Piola 旋度映射)
         */
    void curl(const int cell, const double* lambda, double* values) const;
private:
    /**
         * 参考空间→物理空间的映射信息，basis() 与 curl() 共用。
         */
    struct PrismJacobian {
        double J[3][3];       // 雅可比矩阵 ∂x/∂ξ
        double invJ[3][3];    // 逆雅可比矩阵
        double detJ;          // 行列式
    };

    /**
         * 在指定积分点计算 Jacobian 及其逆。
         * @param cell   单元编号
         * @param lambda 参考坐标 (xi, eta, zeta)
         * @param P      输出：6 个节点的物理坐标
         * @param jac    输出：Jacobian、逆矩阵、行列式
         */
    void computeJacobian(const int cell, const double* lambda,
                         double P[6][3], PrismJacobian& jac) const;

    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::EdgeData<NDIM, int> > d_edge_order;
    tbox::Pointer<pdat::NodeData<NDIM, double> > d_node_coord; // 物理坐标

    tbox::Array<int> d_cell_edge_ext;
    tbox::Array<int> d_cell_edge_idx;
    tbox::Array<int> d_edge_node_ext;
    tbox::Array<int> d_edge_node_idx;
    tbox::Array<int> d_cell_node_ext;
    tbox::Array<int> d_cell_node_idx;
};
}
}
#endif
