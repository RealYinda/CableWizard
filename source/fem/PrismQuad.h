//
// 文件名:      PrismQuad.h
// 描述  :      三棱柱 (Prism) 单元的积分与四次规则管理类，结构参照 TetQuad.h
//

#ifndef included_appu_PrismQuad
#define included_appu_PrismQuad

#ifndef NDIM
#define NDIM 3
#endif

#include <algorithm>
#include <strings.h>
#include <math.h>
#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "CellData.h"
#include "NodeData.h"
#include "GridInfo.h"
namespace JAUMIN {
namespace appu {

/**
 * 类 PrismQuad 负责管理三棱柱单元的数值积分。
 * 包含体积分、面积分（区分三角面和四边形面）以及边积分。
 */
class PrismQuad
{
public:
    /**
     * 结构体 Quad 定义了具体的积分规则数据。
     */
    struct Quad {
        const char *name;
        int dim;
        int order;
        int npoints;
        double *points;
        double *weights;
        int id;
    };

    PrismQuad(const hier::Patch<NDIM>& patch,
              tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
              tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

private:
    /**
     * 验证积分阶数是否受支持。
     */
    void check(int order) const;

    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_volume;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;

    // 拓扑辅助数组
    tbox::Array<int> can_ext, can_idx;
    tbox::Array<int> fan_ext, fan_idx;
    tbox::Array<int> caf_ext, caf_idx;

    // 积分表指针
    Quad** d_quad_table;            // 3D 体积分表
    Quad** d_tri_face_quad_table;   // 三角形面积分表
    Quad** d_quad_face_quad_table;  // 四边形面积分表
    Quad** d_edge_quad_table;       // 1D 边积分表

};
}
}
#endif
