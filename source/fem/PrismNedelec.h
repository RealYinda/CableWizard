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
