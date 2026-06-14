//
// 文件名:      Nedelec.h
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

#ifndef included_appu_triNedelec
#define included_appu_triNedelec

//#if NDIM == 2
//#error FOR 3D CASE ONLY
//#endif

#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"

namespace JAUMIN {
namespace appu {

class triNedelec {
public:
  /**
   * Constructor.
   */
  triNedelec(const hier::Patch<NDIM>& patch,
          tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
          tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

  /**
   * Number of basis functions in a cell.
   */
  int nbas() const;

  /**
   * Dimension of a basis function.
   */
  int dim() const;

  /**
   * Basis function values.
   */
  void basis(const int face, const double* lambda, double* values,const tbox::Array<tbox::Array<double> > &nodeGradBas,int faceorder) const;

  /**
   * Curl of basis functions.
   */
  void curl(const int cell, const double* lambda, double* values,const tbox::Array<tbox::Array<double> > &nodeGradBas,int faceorder) const;

private:
  const hier::Patch<NDIM>& d_patch;
  tbox::Pointer<pdat::EdgeData<NDIM, bool> > d_edge_order;
  tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;
  tbox::Pointer<pdat::NodeData<NDIM, double> > d_node_coord;
  tbox::Array<int> d_face_edge_ext;
  tbox::Array<int> d_face_edge_idx;
  tbox::Array<int> d_edge_node_ext;
  tbox::Array<int> d_edge_node_idx;
  tbox::Array<int> d_face_node_ext;
  tbox::Array<int> d_face_node_idx;
};

}  // namespace appu
}  // namespace JAUMIN

#endif  // included_appu_Nedelec
