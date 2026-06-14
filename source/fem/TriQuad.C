//
// 文件名:      TetQuad.C
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

/* NOTE:
 *   Some of the following codes are borrowed from PHG(Parallel Hierarchical Grid)
 *   lsec.cc.ac.cn/phg */

#include "TriQuad.h"

namespace JAUMIN {
namespace appu {
namespace {

/* 1D */
#define Perm2(a)	a,a
#define Dup2(w)		w
#define Perm11(a)	a,1.-(a),  1.-(a),a
#define Dup11(w)	w,w

/* 2D */
#define Perm3(a)	a,a,a
#define Dup3(w)		w
#define Perm21(a)	a,a,1.-(a)-(a), a,1.-(a)-(a),a, 1.-(a)-(a),a,a
#define Dup21(w)	w,w,w
#define Perm111(a,b)	a,b,1.-(a)-(b), a,1.-(a)-(b),b, \
			b,a,1.-(a)-(b), b,1.-(a)-(b),a, \
			1.-(a)-(b),a,b, 1.-(a)-(b),b,a
#define Dup111(w)	w,w,w,w,w,w

#define Length(wts)	(sizeof(wts) / (sizeof(wts[0])))
#define FLOAT double
#define QUAD TriQuad::Quad

/*--------------------------- 1D cubature rules ------------------------*/
static FLOAT QUAD_1D_P1_wts[] = {
    Dup2(1.L)
};
static FLOAT QUAD_1D_P1_pts[Length(QUAD_1D_P1_wts) * 2] = {
    Perm2(.5L)
};
QUAD QUAD_1D_P1_ = {
    "1D P1",			/* name */
    1,				/* dim */
    1,				/* order */
    Length(QUAD_1D_P1_wts),	/* npoints */
    QUAD_1D_P1_pts,		/* points */
    QUAD_1D_P1_wts,		/* weights */
    -1				/* id */
};

static FLOAT QUAD_1D_P3_wts[] = {
    Dup11(.5L)
};
static FLOAT QUAD_1D_P3_pts[Length(QUAD_1D_P3_wts) * 2] = {
    /* (3 - sqrt(3)) / 6, (3 + sqrt(3)) / 6 */
    Perm11(.2113248654051871177454256097490212L)
};
QUAD QUAD_1D_P3_ = {
    "1D P3",			/* name */
    1,				/* dim */
    3,				/* order */
    Length(QUAD_1D_P3_wts),	/* npoints */
    QUAD_1D_P3_pts,		/* points */
    QUAD_1D_P3_wts,		/* weights */
    -1				/* id */
};

/* 1D quad table */
TriQuad::Quad* edge_quad_table[] = {
    &QUAD_1D_P1_, /* p0 and p1 are the same */
    &QUAD_1D_P1_,
    &QUAD_1D_P1_, /* p2 and p1 are the same */
    &QUAD_1D_P3_
};


/*--------------------------- 2D cubature rules ------------------------*/

static FLOAT QUAD_2D_P1_wts[] = {
    Dup3(1.000000000000000000000000000000000L)
};
static FLOAT QUAD_2D_P1_pts[Length(QUAD_2D_P1_wts) * 3] = {
    /* 1/3 */
    Perm3(.3333333333333333333333333333333333L)
};
QUAD QUAD_2D_P1_ = {
    "2D P1",			/* name */
    2,				/* dim */
    1,				/* order */
    Length(QUAD_2D_P1_wts),	/* npoints */
    QUAD_2D_P1_pts,		/* points */
    QUAD_2D_P1_wts,		/* weights */
    -1				/* id */
};

static FLOAT QUAD_2D_P2_wts[] = {
    /* 1/3 */
    Dup21(.3333333333333333333333333333333333L)
};
static FLOAT QUAD_2D_P2_pts[Length(QUAD_2D_P2_wts) * 3] = {
    /* 1/6 */
    Perm21(.1666666666666666666666666666666667L)
};

//.1666666666666666666666666666666667L,.1666666666666666666666666666666667L,1.-(.1666666666666666666666666666666667L)-(.1666666666666666666666666666666667L), .1666666666666666666666666666666667L,1.-(.1666666666666666666666666666666667L)-(.1666666666666666666666666666666667L),.1666666666666666666666666666666667L, 1.-(.1666666666666666666666666666666667L)-(.1666666666666666666666666666666667L),.1666666666666666666666666666666667L,.1666666666666666666666666666666667L
QUAD QUAD_2D_P2_ = {
    "2D P2",			/* name */
    2,				/* dim */
    2,				/* order */
    Length(QUAD_2D_P2_wts),	/* npoints = 3 */
    QUAD_2D_P2_pts,		/* points */
    QUAD_2D_P2_wts,		/* weights */
    -1				/* id */
};

#if 0
/* Note: this rule has points on the edges */
static FLOAT QUAD_2D_P3_wts[] = {
    /* 1/30 */
    Dup21(.0333333333333333333333333333333333L),
    /* 3/10 */
    Dup21(.3000000000000000000000000000000000L)
};
static FLOAT QUAD_2D_P3_pts[Length(QUAD_2D_P3_wts) * 3] = {
    Perm21(.5000000000000000000000000000000000L),
    /* 1/6 */
    Perm21(.1666666666666666666666666666666667L)
};
#else
static FLOAT QUAD_2D_P3_wts[] = {
    Dup21(.2811498024409796482535143227020770L),
    Dup21(.0521835308923536850798190106312564L)
};
static FLOAT QUAD_2D_P3_pts[Length(QUAD_2D_P3_wts) * 3] = {
    Perm21(.1628828503958919109001618041849063L),
    Perm21(.4779198835675637000000000000000000L)
};
#endif
QUAD QUAD_2D_P3_ = {
    "2D P3",			/* name */
    2,				/* dim */
    3,				/* order */
    Length(QUAD_2D_P3_wts),	/* npoints = 6 */
    QUAD_2D_P3_pts,		/* points */
    QUAD_2D_P3_wts,		/* weights */
    -1				/* id */
};
/* 2D quad table */
TriQuad::Quad* face_quad_table[] = {
    &QUAD_2D_P1_, /* p0 and p1 are the same */
    &QUAD_2D_P1_,
    &QUAD_2D_P2_,
    &QUAD_2D_P3_
};
#undef Perm3
#undef Dup3
#undef Perm21
#undef Dup21
#undef Perm111
#undef Dup111
#undef FLOAT
#undef QUAD
#undef Length


#undef Perm4
#undef Dup4
#undef Perm31
#undef Dup31

} // anoynomous namespace
} // namespace appu
} // namespace JAUMIN


namespace JAUMIN {
namespace appu {

TriQuad::TriQuad(
	const hier::Patch<NDIM>& patch,
	tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
	tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian) :
    d_patch(patch),
    d_cell_volume(cell_volume),
    d_cell_jacobian(cell_jacobian)
{
    d_face_quad_table = &(face_quad_table[0]);
    d_edge_quad_table = &(edge_quad_table[0]);

    patch.getPatchTopology()->getCellAdjacencyNodes(can_ext, can_idx);
    patch.getPatchTopology()->getCellAdjacencyFaces(caf_ext, caf_idx);
    patch.getPatchTopology()->getFaceAdjacencyNodes(fan_ext, fan_idx);
    patch.getPatchTopology()->getCellAdjacencyEdges(cae_ext, cae_idx);
     patch.getPatchTopology()->getEdgeAdjacencyNodes(ean_ext, ean_idx);
}
} // namespace appu
} // namespace JAUMIN
