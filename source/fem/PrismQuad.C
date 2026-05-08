///
/// 文件名:      PrismQuad.C
/// 软件包:
/// 描述  :      三棱柱 (Prism) 单元的积分点信息，格式参照 TetQuad.C
///

#include "PrismQuad.h"
namespace JAUMIN{
namespace appu{
namespace{
// 辅助宏定义，参考 TetQuad.C，用于简化数组的书写
#define Dup6(w) w,w,w,w,w,w
/* ========================================================================
 * 三棱柱 3D 体积分规则 (基于 PrismQuadratureInfo.C)
 * ======================================================================== */

/* * 3D Prism order 1 (代数精度为 1)
 * 1 点积分 (中心点)
 * 坐标: (1/3, 1/3, 0)
 * 权重: 1.0
 */
double quad_prism_3d_p1_points[] = {
  .3333333333333333333333333333333333L, // x
  .3333333333333333333333333333333333L, // y
  0.0000000000000000000000000000000000L  // z
};
double quad_prism_3d_p1_weights[] = {
  1.0000000000000000000000000000000000L
};
PrismQuad::Quad quad_prism_3d_p1 = {
    "Prism-3D-P1",          /* name */
    3,                      /* dim */
    1,                      /* order */
    1,                      /* npoints */
    quad_prism_3d_p1_points, /* points */
    quad_prism_3d_p1_weights, /* weights */
    -1
};
/* * 3D Prism order 2 (代数精度为 2 及以上)
 * 6 点积分
 * 权重: 均为 1/6
 * Z坐标: +/- 1/sqrt(3) ≈ +/- 0.57735026918962576451
 */
double quad_prism_3d_p2_points[] = {
    // 下半层积分点 (zeta = -1 / sqrt(3))
    .1666666666666666666666666666666667L, // 1/6
    .1666666666666666666666666666666667L, // 1/6
   -.5773502691896257645091487805019574L, // -1/sqrt(3)

    .6666666666666666666666666666666667L, // 2/3
    .1666666666666666666666666666666667L, // 1/6
   -.5773502691896257645091487805019574L,

    .1666666666666666666666666666666667L, // 1/6
    .6666666666666666666666666666666667L, // 2/3
   -.5773502691896257645091487805019574L,

    // 上半层积分点 (zeta = +1 / sqrt(3))
    .1666666666666666666666666666666667L,
    .1666666666666666666666666666666667L,
    .5773502691896257645091487805019574L, // +1/sqrt(3)

    .6666666666666666666666666666666667L,
    .1666666666666666666666666666666667L,
    .5773502691896257645091487805019574L,

    .1666666666666666666666666666666667L,
    .6666666666666666666666666666666667L,
    .5773502691896257645091487805019574L
};

double quad_prism_3d_p2_weights[] = {
    Dup6(.1666666666666666666666666666666667L) // 1/6
};

PrismQuad::Quad quad_prism_3d_p2 = {
    "Prism-3D-P2",          /* name */
    3,                      /* dim */
    2,                      /* order */
    6,                      /* npoints */
    quad_prism_3d_p2_points, /* points */
    quad_prism_3d_p2_weights, /* weights */
    -1
};

/* 三棱柱体积分表 */
PrismQuad::Quad* quad_table[] = {
    &quad_prism_3d_p1, /* P0 和 P1 共用 1阶积分规则 */
    &quad_prism_3d_p1,
    &quad_prism_3d_p2, /* P2 及以上使用 6点积分规则 */
    &quad_prism_3d_p2  /* 假设更高阶暂时沿用这个，需要时可扩展 */
};


/* ========================================================================
 * TODO: 如果需要为三棱柱的“面(四边形面/三角形面)”和“边”做积分
 * 请参考 TetQuad.C 在这里补充 face_quad_table 和 edge_quad_table。
 * 目前暂留空表以防止编译报错。
 * ======================================================================== */
// 辅助宏定义
#define Perm2(a)    a,a
#define Dup2(w)     w
#define Perm11(a)   a,1.-(a),  1.-(a),a
#define Dup11(w)    w,w

#define Perm3(a)    a,a,a
#define Dup3(w)     w
#define Perm21(a)   a,a,1.-(a)-(a), a,1.-(a)-(a),a, 1.-(a)-(a),a,a
#define Dup21(w)    w,w,w

#define Length(wts) (sizeof(wts) / (sizeof(wts[0])))
#define FLOAT double

static FLOAT QUAD_1D_P1_wts[] = { Dup2(1.L) };
static FLOAT QUAD_1D_P1_pts[Length(QUAD_1D_P1_wts) * 2] = { Perm2(.5L) };
PrismQuad::Quad QUAD_1D_P1_ = {
  "1D-P1", 1, 1, Length(QUAD_1D_P1_wts), QUAD_1D_P1_pts, QUAD_1D_P1_wts, -1
};

static FLOAT QUAD_1D_P3_wts[] = { Dup11(.5L) };
static FLOAT QUAD_1D_P3_pts[Length(QUAD_1D_P3_wts) * 2] = {
    Perm11(.2113248654051871177454256097490212L) // (3 - sqrt(3)) / 6, (3 + sqrt(3)) / 6
};
PrismQuad::Quad QUAD_1D_P3_ = {
  "1D-P3", 1, 3, Length(QUAD_1D_P3_wts), QUAD_1D_P3_pts, QUAD_1D_P3_wts, -1
};

/* 边积分表 */
PrismQuad::Quad* edge_quad_table[] = {
  &QUAD_1D_P1_, /* P0 and P1 */
  &QUAD_1D_P1_,
  &QUAD_1D_P1_, /* P2 */
  &QUAD_1D_P3_  /* P3 */
};

/* ========================================================================
 * 2D Triangle Face 积分规则 (三角形面积分，基于面重心坐标)
 * ======================================================================== */
static FLOAT QUAD_TRI_P1_wts[] = { Dup3(1.000000000000000000000000000000000L) };
static FLOAT QUAD_TRI_P1_pts[Length(QUAD_TRI_P1_wts) * 3] = { Perm3(.3333333333333333333333333333333333L) };
PrismQuad::Quad QUAD_TRI_P1_ = {
    "2D-TRI-P1", 2, 1, Length(QUAD_TRI_P1_wts), QUAD_TRI_P1_pts, QUAD_TRI_P1_wts, -1
};

static FLOAT QUAD_TRI_P2_wts[] = { Dup21(.3333333333333333333333333333333333L) };
static FLOAT QUAD_TRI_P2_pts[Length(QUAD_TRI_P2_wts) * 3] = { Perm21(.1666666666666666666666666666666667L) };
PrismQuad::Quad QUAD_TRI_P2_ = {
    "2D-TRI-P2", 2, 2, Length(QUAD_TRI_P2_wts), QUAD_TRI_P2_pts, QUAD_TRI_P2_wts, -1
};

/* 三角形面积分表 */
PrismQuad::Quad* tri_face_quad_table[] = {
    &QUAD_TRI_P1_,
    &QUAD_TRI_P1_,
    &QUAD_TRI_P2_,
    &QUAD_TRI_P2_ // 若有高阶可以继续补充
};


/* ========================================================================
 * 2D Quadrilateral Face 积分规则 (四边形面积分，基于自然坐标 [-1,1]x[-1,1] 或 [0,1]x[0,1])
 * 注意：此处假设使用标准的二维 Gauss 积分。具体点格式需要与你 bary2dTo3d 函数匹配。
 * ======================================================================== */
// 1点积分 (中心点) - 对应1阶精度
static FLOAT QUAD_QUAD_P1_wts[] = { 4.0000000000000000000000000000000000L };
static FLOAT QUAD_QUAD_P1_pts[] = { 0.0L, 0.0L }; // 假设四边形局部坐标在 [-1, 1]
PrismQuad::Quad QUAD_QUAD_P1_ = {
    "2D-QUAD-P1", 2, 1, 1, QUAD_QUAD_P1_pts, QUAD_QUAD_P1_wts, -1
};

// 4点积分 (2x2 Gauss点) - 对应2阶精度以上
static FLOAT QUAD_QUAD_P2_wts[] = {
    1.0L, 1.0L, 1.0L, 1.0L
};
static FLOAT QUAD_QUAD_P2_pts[] = {
    -.57735026918962576451L, -.57735026918962576451L, // (-1/sqrt(3), -1/sqrt(3))
     .57735026918962576451L, -.57735026918962576451L, // (+1/sqrt(3), -1/sqrt(3))
     .57735026918962576451L,  .57735026918962576451L, // (+1/sqrt(3), +1/sqrt(3))
    -.57735026918962576451L,  .57735026918962576451L  // (-1/sqrt(3), +1/sqrt(3))
};
PrismQuad::Quad QUAD_QUAD_P2_ = {
    "2D-QUAD-P2", 2, 2, 4, QUAD_QUAD_P2_pts, QUAD_QUAD_P2_wts, -1
};

/* 四边形面积分表 */
PrismQuad::Quad* quad_face_quad_table[] = {
    &QUAD_QUAD_P1_,
    &QUAD_QUAD_P1_,
    &QUAD_QUAD_P2_,
    &QUAD_QUAD_P2_
};

#undef Perm2
#undef Dup2
#undef Perm11
#undef Dup11
#undef Perm3
#undef Dup3
#undef Perm21
#undef Dup21
#undef Length
#undef FLOAT

#undef Dup6

} // anonymous namespace
} // namespace appu
} // namespace JAUMIN


namespace JAUMIN {
namespace appu {

// 构造函数：获取 patch 拓扑结构以及单元 Jacobian 映射等信息
PrismQuad::PrismQuad(
    const hier::Patch<NDIM>& patch,
    tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
    tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian) :
  d_patch(patch),
  d_cell_volume(cell_volume),
  d_cell_jacobian(cell_jacobian)
{
  // 将指针挂载到匿名空间的静态数组上
  d_quad_table = &(quad_table[0]);
  d_tri_face_quad_table = &(tri_face_quad_table[0]);
  d_quad_face_quad_table = &(quad_face_quad_table[0]);
  d_edge_quad_table = &(edge_quad_table[0]);

  patch.getPatchTopology()->getCellAdjacencyNodes(can_ext, can_idx);
  patch.getPatchTopology()->getCellAdjacencyFaces(caf_ext, caf_idx);
  patch.getPatchTopology()->getFaceAdjacencyNodes(fan_ext, fan_idx);
}

} // namespace appu
} // namespace JAUMIN
