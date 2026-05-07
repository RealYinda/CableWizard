//
// 文件名:      TetQuad.h
// 软件包:      
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :      
//

#ifndef included_appu_TetQuad
#define included_appu_TetQuad

#if NDIM == 2
#error FOR 3D CASE ONLY
#endif

#include <algorithm>
#include <strings.h>
#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"
#include "GridInfo.h"

namespace JAUMIN {
namespace appu {

class TetQuad
{
public:
    TetQuad(const hier::Patch<NDIM>& patch, 
            tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
	    tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

    template<class ShapeFunc>
    void quadBasDotBas(const int cell, 
	               const ShapeFunc* shapefunc, 
	               const int order, 
		       double *value) const
    {
	check(order);
	int nbas = shapefunc->nbas();
	int dim = shapefunc->dim();
        double (*bb)[nbas] = (double (*)[nbas])value;
	const Quad *quad = d_quad_table[order];

	double vol = (*d_cell_volume)(0, cell);
	bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));

	for(int n = 0; n < quad->npoints; n++) {
	    double phi[nbas][dim];
	    shapefunc->basis(cell,
		             quad->points + n * (NDIM + 1),
			     &(phi[0][0]));
	    for(int i = 0; i < nbas; i++) {
		for(int j = 0; j < nbas; j++) {
		    double v = 0.;
		    for(int d = 0; d < dim; d++)
			v += phi[i][d] * phi[j][d];
		    v *= quad->weights[n];
		    bb[i][j] += v;
		}
	    }
	}
	for(int i = 0; i < nbas; i++)
	    for(int j = 0; j < nbas; j++)
		bb[i][j] *= vol;
    }


    /**
     * 将2d重心坐标lambda_2d映射到第cell个单元的第face个面上。
     */
    void bary2dTo3d(const double* lambda_2d, int cell, int face, double* lambda_3d) const {
        int patch_face = caf_idx[caf_ext[cell] + face];
        int cnt = 0;
        for (int i=0; i<4; i++) {
            lambda_3d[i] = 0;
        }
        // 对面的顶点进行循环
        for (int j = fan_ext[patch_face]; j<fan_ext[patch_face+1]; j++) {
            int face_node = fan_idx[j];
            // 对单元的结点进行循环，得到face_node的局部编号
            for (int i=can_ext[cell]; i<can_ext[cell+1]; i++) {
                int cell_node = can_idx[i];
                int local_node = i-can_ext[cell]; // local_node为面的局部编号
                if (cell_node == face_node) {
                    lambda_3d[local_node] = lambda_2d[cnt++];
                }
            }
        }
    }


    template<class ShapeFunc>
    void faceQuadBasDotBas(const int cell, 
                           const int face,
	               const ShapeFunc* shapefunc, 
	               const int order, 
		       double *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        double (*bb)[nbas] = (double (*)[nbas])value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);


        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw)) / 2.0;


        for(int n = 0; n < quad->npoints; n++) {
            // transform 2d barycentic to 3d
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n * (NDIM + 1 - 1)), cell, face, bary3d);

            // get the shape functions' value
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));

            // compute the stiff matrix of face 'face' of element 'cell'
            for(int i = 0; i < nbas; i++) {
                for(int j = 0; j < nbas; j++) {
                    double bas_i[NDIM], bas_j[NDIM];
                    crossProduct(NDIM, outer_normal, phi[i], bas_i);
                    crossProduct(NDIM, outer_normal, phi[j], bas_j);
                    double v = dotProduct(NDIM, bas_i, bas_j);
                    v *= quad->weights[n];
                    bb[i][j] += v;
                }
            }
        }
        for(int i = 0; i < nbas; i++)
            for(int j = 0; j < nbas; j++)
                bb[i][j] *= area;
    }

    template<class ShapeFunc, class TYPE, class Function>
    void faceQuadFunctionDotBas(const int cell,
                                const int face,
                                const Function& func,
                                const ShapeFunc* shapefunc,
                                const int order,
                                TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        TYPE *bb = value;
        const Quad *quad = d_face_quad_table[order];
        assert(face < NDIM+1 && face >= 0);

        //double vol = (*d_cell_volume)(0, cell);
        bzero(&(bb[0]), nbas * sizeof(*value));

        double outer_normal[NDIM], outer_normal_raw[NDIM];
        // get the id of face in this patch.
        int patch_face = caf_idx[caf_ext[cell] + face];
        outerNormal(
                can_ext[cell+1]-can_ext[cell],
                can_idx.getPointer() + can_ext[cell],
                fan_ext[patch_face+1]-fan_ext[patch_face],
                fan_idx.getPointer() + fan_ext[patch_face],
                d_patch.getNumberOfNodes(1),
                d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                NDIM,
                outer_normal,
                outer_normal_raw);
        double area = sqrt(dotProduct(NDIM, outer_normal_raw, outer_normal_raw)) / 2.0;


        for(int n = 0; n < quad->npoints; n++) {
            // transform this 2d barycentric quadrature point to 3d barycentic on face 'face'
            double bary3d[NDIM+1];
            bary2dTo3d(quad->points + (n * (NDIM + 1 - 1)), cell, face, bary3d);

            // get the space coordinate of this quadrature point
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               bary3d,
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);

            // get the 'nbas' shape functions' value at this quadrature point
            double phi[nbas][dim];
            shapefunc->basis(cell,
                             bary3d,
                             &(phi[0][0]));
            // compute the load of this element
            for(int i = 0; i < nbas; i++) {
                TYPE v = dotProduct(NDIM, funcvalues, phi[i]);
                v *= quad->weights[n];
                bb[i] += v;
            }
        }
        for(int i = 0; i < nbas; i++)
                bb[i] *= area;
    }


    template<class ShapeFunc>
    void quadGradBasDotGradBas(const int cell, 
	                       const ShapeFunc* shapefunc, 
	                       const int order, 
		               double *value) const
    {
	check(order);
	TBOX_ASSERT(shapefunc->dim() == 1); 

	int nbas = shapefunc->nbas();
        double (*gg)[nbas] = (double (*)[nbas])value;
	const Quad *quad = d_quad_table[order];

	double vol = (*d_cell_volume)(0, cell);
	bzero(&(gg[0][0]), nbas * nbas * sizeof(*value));

	for(int n = 0; n < quad->npoints; n++) {
	    double gradphi[nbas][NDIM];
	    shapefunc->gradient(cell,
		                quad->points + n * (NDIM + 1),
			        &(gradphi[0][0]));
	    for(int i = 0; i < nbas; i++) {
		for(int j = 0; j < nbas; j++) {
		    double v = 0.;
		    for(int d = 0; d < NDIM; d++)
			v += gradphi[i][d] * gradphi[j][d];
		    v *= quad->weights[n];
		    gg[i][j] += v;
		}
	    }
	}
	for(int i = 0; i < nbas; i++)
	    for(int j = 0; j < nbas; j++)
		gg[i][j] *= vol;
    }

    template<class ShapeFunc>
    void quadCurlBasDotCurlBas(const int cell, 
	                       const ShapeFunc* shapefunc, 
	                       const int order, 
		               double *value) const
    {
	check(order);
	TBOX_ASSERT(shapefunc->dim() == NDIM); 

	int nbas = shapefunc->nbas();
        double (*cc)[nbas] = (double (*)[nbas])value;
	const Quad *quad = d_quad_table[order];

	double vol = (*d_cell_volume)(0, cell);
	bzero(&(cc[0][0]), nbas * nbas * sizeof(*value));

	for(int n = 0; n < quad->npoints; n++) {
	    double curlphi[nbas][NDIM];
	    shapefunc->curl(cell,
		            quad->points + n * (NDIM + 1),
			    &(curlphi[0][0]));
	    for(int i = 0; i < nbas; i++) {
		for(int j = 0; j < nbas; j++) {
		    double v = 0.;
		    for(int d = 0; d < NDIM; d++)
			v += curlphi[i][d] * curlphi[j][d];
		    v *= quad->weights[n];
		    cc[i][j] += v;
		}
	    }
	}
	for(int i = 0; i < nbas; i++)
	    for(int j = 0; j < nbas; j++)
		cc[i][j] *= vol;
    }

    template<class ShapeFunc, class TYPE>
    void quadDofDotBas(const int cell, 
	               const ShapeFunc* shapefunc, 
	               const int order, 
               const TYPE *dof_value,
               TYPE *value) const
    {
	check(order);
	int nbas = shapefunc->nbas();
	int dim = shapefunc->dim();
	const Quad *quad = d_quad_table[order];

	double vol = (*d_cell_volume)(0, cell);
	bzero(value, nbas * sizeof(*value));

	for(int n = 0; n < quad->npoints; n++) {
	    double phi[nbas][dim];
	    shapefunc->basis(cell,
		             quad->points + n * (NDIM + 1),
			     &(phi[0][0]));
        for(int i = 0; i < nbas; i++) {
            TYPE v = 0.;
            for(int j = 0; j < nbas; j++) {
                TYPE u = 0.;
                for(int d = 0; d < dim; d++)
                    u += phi[i][d] * phi[j][d];
                v += u * dof_value[j];
            }
            value[i] += v * quad->weights[n];
	    }
	}
	for(int i = 0; i < nbas; i++)
	  value[i] *= vol;
    }


    template<class ShapeFunc, class TYPE, class Function>
    void quadFunctionDotBas(
            const int cell,
            const Function& func,
            const ShapeFunc* shapefunc,
            const int order,
            TYPE *value) const
    {
        check(order);
        int nbas = shapefunc->nbas();
        int dim = shapefunc->dim();
        const Quad *quad = d_quad_table[order];

        double vol = (*d_cell_volume)(0, cell);
        bzero(value, nbas * sizeof(*value));

        for(int n = 0; n < quad->npoints; n++) {
            double phi[nbas][dim];
            double space_coords[NDIM];
            barycentricToSpace(NDIM,
                               d_patch.getNumberOfNodes(1),
                               d_patch.getPatchGeometry()->getNodeCoordinates()->getPointer(),
                               can_idx.getPointer() + can_ext[cell],
                               quad->points + n *(NDIM + 1),
                               space_coords);
            TYPE funcvalues[NDIM];
            func(space_coords[0],space_coords[1],space_coords[2], funcvalues);
            shapefunc->basis(cell,
                             quad->points + n * (NDIM + 1),
                             &(phi[0][0]));
            for(int i = 0; i < nbas; i++) {
                TYPE prod = dotProduct(NDIM, funcvalues, phi[i]);
                value[i] += prod * quad->weights[n];
            }
        }
        for(int i = 0; i < nbas; i++)
            value[i] *= vol;
    }

    template<class ShapeFunc, class TYPE>
    void quadCurlDofDotCurlDof(const int cell, 
	                       const ShapeFunc* shapefunc, 
	                       const int order, 
                       const TYPE *dof_u,
                       const TYPE *dof_v,
                       TYPE *value) const
    {
	check(order);
	int nbas = shapefunc->nbas();
	int dim = shapefunc->dim();
	const Quad *quad = d_quad_table[order];

	double vol = (*d_cell_volume)(0, cell);
        *value = 0.;

	for(int n = 0; n < quad->npoints; n++) {
	    double curlphi[nbas][NDIM];
	    shapefunc->curl(cell,
		            quad->points + n * (NDIM + 1),
			    &(curlphi[0][0]));

            TYPE uv[NDIM], vv[NDIM];
            for(int d = 0; d < NDIM; d++)
                uv[d] = vv[d] = 0.;

	    for(int i = 0; i < nbas; i++) {
                for(int d = 0; d < dim; d++) {
                    uv[d] += dof_u[i] * curlphi[i][d];
                    vv[d] += dof_v[i] * curlphi[i][d];
                }
            }

            TYPE q = dotProduct(dim, uv, vv);

            //TYPE q = 0.;
            //for(int d = 0; d < dim; d++)
            //    q += uv[d] * vv[d];

            *value += q * quad->weights[n];
	}

        *value *= vol;
    }

public:
    struct Quad {
	const char *name;
	int dim;
	int order;
	int npoints;
	double *points;
	double *weights;
    int id;  
    };

private:
    void check(int order) const {
	const int max_order = 3;
	if(order > max_order)
	    TBOX_ERROR("Quadrature rules of order " 
		       << order << " for tetrahedra not implememted.\n");
    }

    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_volume;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;

    tbox::Array<int> can_ext, can_idx;
    tbox::Array<int> fan_ext, fan_idx;
    tbox::Array<int> caf_ext, caf_idx;

    Quad** d_quad_table;
    Quad** d_face_quad_table;
    Quad** d_edge_quad_table;
};

} // namespace appu
} // namespace JAUMIN

#endif // included_appu_TetQuad
