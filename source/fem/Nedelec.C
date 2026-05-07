//
// 文件名:      Nedelec.C
// 软件包:      
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :      
//

#include "Nedelec.h"

namespace JAUMIN {
namespace appu {

Nedelec::Nedelec(const hier::Patch<NDIM>& patch, 
                 tbox::Pointer<pdat::EdgeData<NDIM, int> > edge_order,
                 tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian) :
    d_patch(patch),
    d_edge_order(edge_order),
    d_cell_jacobian(cell_jacobian)
{
    patch.getPatchTopology()->getCellAdjacencyEdges(d_cell_edge_ext, d_cell_edge_idx);
    patch.getPatchTopology()->getEdgeAdjacencyNodes(d_edge_node_ext, d_edge_node_idx);
    patch.getPatchTopology()->getCellAdjacencyNodes(d_cell_node_ext, d_cell_node_idx);
}

int Nedelec::nbas() const
{ return NDIM == 2 ? 3 : 6; }

int Nedelec::dim() const
{ return NDIM; }

void Nedelec::basis(const int cell, const double* lambda, double* values) const
/*
 * \phi_ij =  \lambda_i \grad \lambda_j - \lambda_j \grad \lambda_i
 */
{
    double (*nabla)[NDIM + 1] = (double (*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
    double (*phi)[NDIM] = (double (*)[NDIM])values;

    int nedge = NDIM == 2 ? 3 : 6;
    for(int ie = 0; ie < nedge; ie++) {
        int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
	int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
	int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];
	int in_order = (*d_edge_order)(0, edge);
	if(!in_order) {
	    int t = v0;
	    v0 = v1; v1 = t;
	}

        int lv0, lv1;
        lv0 = lv1 = -1;
        for(int n = 0; n < NDIM + 1; n++) {
            int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
            if(nidx == v0) lv0 = n;
            if(nidx == v1) lv1 = n;
        }

	phi[ie][0] = lambda[lv0] * nabla[lv1][0] - lambda[lv1] * nabla[lv0][0];
	phi[ie][1] = lambda[lv0] * nabla[lv1][1] - lambda[lv1] * nabla[lv0][1];
	phi[ie][2] = lambda[lv0] * nabla[lv1][2] - lambda[lv1] * nabla[lv0][2];
    }
}

void Nedelec::curl(const int cell, const double* lambda, double* values) const
/*
 * \curl \phi_ij =  2\grad\lambda_i x \grad\lambda_j
 */
{
    NULL_USE(lambda);
    double (*nabla)[NDIM + 1] = (double (*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
    double (*curlphi)[NDIM] = (double (*)[NDIM])values;
    int nedge = NDIM == 2 ? 3 : 6;
    for(int ie = 0; ie < nedge; ie++) {
        int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
	int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
	int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];
	int in_order = (*d_edge_order)(0, edge);
	if(!in_order) {
	    int t = v0;
	    v0 = v1; v1 = t;
	}

        int lv0, lv1;
        lv0 = lv1 = -1;
        for(int n = 0; n < NDIM + 1; n++) {
            int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
            if(nidx == v0) lv0 = n;
            if(nidx == v1) lv1 = n;
        }

	double a = nabla[lv0][0];
	double b = nabla[lv0][1];
	double c = nabla[lv0][2];
	double d = nabla[lv1][0];
	double e = nabla[lv1][1];
	double f = nabla[lv1][2];

	curlphi[ie][0] = (b * f - e * c) * 2.;
	curlphi[ie][1] = (c * d - f * a) * 2.;
	curlphi[ie][2] = (a * e - d * b) * 2.;
    }
}

} // namespace appu
} // namespace JAUMIN

