//
// 文件名:      TetGeom.h
// 软件包:      
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :      
//

#ifndef included_appu_TetGeom
#define included_appu_TetGeom

#if NDIM == 2
#error FOR 3D CASE ONLY
#endif

#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "NodeData.h"

namespace JAUMIN {
namespace appu {

class TetGeom
{
public:
    /**
     * Constructor.
     */
    TetGeom(const hier::Patch<NDIM>& patch);

    /**
     * Computes the Jabobian
     *  J = D\lambda / Dx = [ D\lambda_0/Dx, D\lambda_0/Dy, D\lambda_0/Dz, c0;
     *                        D\lambda_1/Dx, D\lambda_1/Dy, D\lambda_1/Dz, c1;
     *                        D\lambda_2/Dx, D\lambda_2/Dy, D\lambda_2/Dz, c2;
     *                        D\lambda_3/Dx, D\lambda_3/Dy, D\lambda_3/Dz, c3 ]
     * of the element.
     */
    void jacobian(const int cell, double *value) const;

    /**
     * Computes volume of a tetrahedron.
     */
    double volume(const int cell) const;

    /**
     * Computes normal of a triangular face.
     */
    void normal(int face, double *n) const;

    /**
     * Computes area of a triangular face.
     */
    double area(int face) const;

private:
    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::NodeData<NDIM, double> > d_node_coord;
    tbox::Array<int> d_cell_node_ext;
    tbox::Array<int> d_cell_node_idx;
    tbox::Array<int> d_face_node_ext;
    tbox::Array<int> d_face_node_idx;
};

} // namespace appu
} // namespace JAUMIN

#endif // included_appu_TetGeom
