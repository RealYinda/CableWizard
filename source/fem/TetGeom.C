//
// 文件名:      TetGeom.C
// 软件包:      
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :      
//

#include <math.h>
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "NodeData.h"

#include "TetGeom.h"

namespace JAUMIN {
namespace appu {

TetGeom::TetGeom(const hier::Patch<NDIM>& patch) :
    d_patch(patch)
{
    patch.getPatchTopology()->getCellAdjacencyNodes(d_cell_node_ext, d_cell_node_idx);
    patch.getPatchTopology()->getFaceAdjacencyNodes(d_face_node_ext, d_face_node_idx);
    d_node_coord = patch.getPatchGeometry()->getNodeCoordinates();
}

void TetGeom::jacobian(const int cell, double *value) const
/* computes the Jabobian
 *  J = D\lambda / Dx = [ D\lambda_0/Dx, D\lambda_0/Dy, D\lambda_0/Dz, c0;
 *                        D\lambda_1/Dx, D\lambda_1/Dy, D\lambda_1/Dz, c1;
 *                        D\lambda_2/Dx, D\lambda_2/Dy, D\lambda_2/Dz, c2;
 *                        D\lambda_3/Dx, D\lambda_3/Dy, D\lambda_3/Dz, c3 ]
 * of the element, by inverting the matrix
 *
 *  [x0, x1, x2, x3; y0, y1, y2, y3; z0, z1, z2, z3; 1 1 1 1].
 *
 * Note: lambda = J[0..NDIM][0..NDIM-1] * x + J[0..NDIM][NDIM]
 * Author: Linbo Zhang(zlb@lsec.cc.ac.cn)
 */
{
    int i, j, k;
    int n = NDIM + 1;
    double d, t, a[NDIM + 1][2 * (NDIM + 1)];
    double (*jacobian)[NDIM + 1] = (double (*)[NDIM + 1])value;
    int *node = &(d_cell_node_idx[d_cell_node_ext[cell]]);

    /* right hand side (I) */
    for (j = 0; j < n; j++) 
	for (i = 0; i < n; i++) 
	    a[j][n + i] = (i == j) ? 1.0 : 0.0;

    /* coefficient matrix */
    for (j = 0; j < n; j++) {
	for (i = 0; i < NDIM; i++)
	    a[j][i] = (*d_node_coord)(i, node[j]);
	a[j][NDIM] = 1.0;
    }

    /* solve linear systems of equations */
    for (j = 0; j < n; j++) {
	/* find pivot */
	d = fabs(a[j][j]);
	k = j;
	for (i = j + 1; i < n; i++)
	    if (d < (t = fabs(a[i][j]))) {
		d = t;
		k = i;
	    }
	if (k != j) /* swap rows j and k */
	    for (i = j; i < n + n; i++) {
		t = a[j][i];
		a[j][i] = a[k][i];
		a[k][i] = t;
	    }
	d = 1.0 / ((t = a[j][j]) == 0. ? 1. : t);
	for (i = j + 1; i < n + n; i++)
		a[j][i] *= d;
	for (i = j + 1; i < n; i++) {
		d = a[i][j];
		for (k = j + 1; k < n + n; k++)
		a[i][k] -= d * a[j][k];
	}
    }

    for (j = n - 2; j >= 0; j--)
	for (i = j + 1; i < n; i++) {
	    d = a[j][i];
	    for (k = n; k < n + n; k++)
	    a[j][k] -= d * a[i][k];
	}

    /* Jacobian = trans(a[0..3, 4..7]) */
    for (j = 0; j < n; j++)
	for (i = 0; i < n; i++)
	    jacobian[j][i] = a[i][j + n];
}

double TetGeom::volume(const int cell) const
/* computes volume of a tetrahedron
 *-----------------------------------------------------------------------
 * General formula for the volume of a k-simplex in n-dimensions:
 *	V = \sqrt(|\det(W * W^T)|)/k!
 * where:
 *	W = [v1-v0,v2-v0,..,vk-v0]
 * For k=n=2:
 *		| x1-x0 y1-y0 |
 *	V = det | x2-x0 y2-y0 | / 2
 * For k=n=3:
 *		| x1-x0 y1-y0 z1-z0 |
 *	V = det | x2-x0 y2-y0 z2-z0 | / 6
 *		| x3-x0 y3-y0 z3-z0 |
 *
 * 3D: V=S*h/3 where S is the area of the base, h the height 
 *----------------------------------------------------------------------
 * Author: Linbo Zhang(zlb@lsec.cc.ac.cn)
 */
{
    double vol;
    double oned6 = 1.0/(double)6.0;
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

    int *node = &(d_cell_node_idx[d_cell_node_ext[cell]]);

    x0 = (*d_node_coord)(0, node[0]);
    y0 = (*d_node_coord)(1, node[0]);
    z0 = (*d_node_coord)(2, node[0]);
    x1 = (*d_node_coord)(0, node[1]) - x0;
    y1 = (*d_node_coord)(1, node[1]) - y0;
    z1 = (*d_node_coord)(2, node[1]) - z0;
    x2 = (*d_node_coord)(0, node[2]) - x0;
    y2 = (*d_node_coord)(1, node[2]) - y0;
    z2 = (*d_node_coord)(2, node[2]) - z0;
    x3 = (*d_node_coord)(0, node[3]) - x0;
    y3 = (*d_node_coord)(1, node[3]) - y0;
    z3 = (*d_node_coord)(2, node[3]) - z0;

    vol = fabs(x1*y2*z3 + x2*y3*z1 + y1*z2*x3 -
		(z1*y2*x3 + y1*x2*z3 + z2*y3*x1)) * oned6;
    return vol;
}

void TetGeom::normal(const int face, double* n) const
/* computes normal of a triangular face
 * Author: Linbo Zhang(zlb@lsec.cc.ac.cn) */
{
    double d, a[3], b[3];
    int v0, v1, v2;
    const double *p0, *p1, *p2;

    /* indices of the three vertices on the face */
    v0 = d_face_node_idx[d_face_node_ext[face] + 0];
    v1 = d_face_node_idx[d_face_node_ext[face] + 1];
    v2 = d_face_node_idx[d_face_node_ext[face] + 2];

    p0 = &((*d_node_coord)(0, v0));
    p1 = &((*d_node_coord)(0, v1));
    p2 = &((*d_node_coord)(0, v2));

    a[0] = p1[0] - p0[0]; a[1] = p1[1] - p0[1]; a[2] = p1[2] - p0[2];
    b[0] = p2[0] - p0[0]; b[1] = p2[1] - p0[1]; b[2] = p2[2] - p0[2];

    n[0] = a[1] * b[2] - a[2] * b[1];
    n[1] = a[2] * b[0] - a[0] * b[2];
    n[2] = a[0] * b[1] - a[1] * b[0];
    d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    d = 1.0 / (d == 0.0 ? 1.0 : sqrt(d));
    n[0] *= d; n[1] *= d; n[2] *= d;
}

double TetGeom::area(const int face) const
/* computes area of a triangular face
 * Author: Linbo Zhang(zlb@lsec.cc.ac.cn) */
{
    double c, d, a[3], b[3], area;
    int v0, v1, v2;
    const double *p0, *p1, *p2;

    /* indices of the three vertices on the face */
    v0 = d_face_node_idx[d_face_node_ext[face] + 0];
    v1 = d_face_node_idx[d_face_node_ext[face] + 1];
    v2 = d_face_node_idx[d_face_node_ext[face] + 2];

    p0 = &((*d_node_coord)(0, v0));
    p1 = &((*d_node_coord)(0, v1));
    p2 = &((*d_node_coord)(0, v2));

    a[0] = p1[0] - p0[0]; a[1] = p1[1] - p0[1]; a[2] = p1[2] - p0[2];
    b[0] = p2[0] - p0[0]; b[1] = p2[1] - p0[1]; b[2] = p2[2] - p0[2];

    /* compute the area (Heron's formula) 
       (http://mathworld.wolfram.com/Circumradius.html) */
    a[0] = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[1] = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    b[0] = p2[0] - p1[0]; b[1] = p2[1] - p1[1]; b[2] = p2[2] - p1[2];
    a[2] = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    d = (a[0] + a[1] + a[2]) * 0.5;

    c = d * (d - a[0]) * (d - a[1]) * (d - a[2]);
    if (c <= 0.) {
	TBOX_WARNING("bad mesh: degenerated face found.\n");
	area = 0.;
    }
    else {
	area = sqrt(c);
    }

    return area;
}

} // namespace appu
} // namespace JAUMIN
