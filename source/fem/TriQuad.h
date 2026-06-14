//
// 文件名:      TetQuad.h
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

#ifndef included_appu_TriQuad
#define included_appu_TriQuad

#include "Array.h"
#include "CellData.h"
#include "EdgeData.h"
#include "GridInfo.h"
#include "NodeData.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchTopology.h"
#include "Pointer.h"
#include <algorithm>
#include <strings.h>

namespace JAUMIN {
namespace appu {
class TriQuad {
public:
  TriQuad(const hier::Patch<NDIM> &patch, tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
          tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

  /****************************************************************************
   * 结点基函数与结点基函数的点积
   ****************************************************************************/
  template <class TYPE>
  void nodeBasdotBas(int order, const bool raw_order, double area, double *value) {
    const Quad *quad = d_face_quad_table[order];
    const int nbas = 3;
    double(*cc)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    int p0 = 0, p1 = 1, p2 = 2;
    if (!raw_order) {
      p1 = 2;
      p2 = 1;
    }
    int perm[3] = {p0, p1, p2};
    // 初始化矩阵
    std::fill(&(cc[0][0]), &(cc[0][0]) + nbas * nbas, 0.0);
    for (int n = 0; n < quad->npoints; n++) {
      const double *p = quad->points + n * 3;
      double w = quad->weights[n] * area;
      // 提取重映射后的形函数值
      double phi[nbas] = {p[perm[0]], p[perm[1]], p[perm[2]]};
      for (int i = 0; i < nbas; i++) {
        double phi_i_w = phi[i] * w;
        for (int j = 0; j < nbas; j++) {
          cc[i][j] += phi_i_w * phi[j];
        }
      }
    }
  }

  /****************************************************************************
   * 结点基函数梯度与结点基函数梯度的点积
   ****************************************************************************/
  template <class TYPE>
  void nodeGraddotGrad(int order, const bool raw_order, double area,
                       const tbox::Array<tbox::Array<double> > &nodeGradBas, double *value) {
    const Quad *quad = d_face_quad_table[order];
    const int nbas = 3;
    double(*cc)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    // 1. 初始化矩阵
    std::fill(&(cc[0][0]), &(cc[0][0]) + nbas * nbas, 0.0);
    int p[3] = {0, 1, 2};
    if (!raw_order) {
      p[1] = 2;
      p[2] = 1;
    }
    double grad[nbas][2];
    for (int i = 0; i < nbas; ++i) {
      grad[i][0] = nodeGradBas[p[i]][0];
      grad[i][1] = nodeGradBas[p[i]][1];
    }
    for (int n = 0; n < quad->npoints; n++) {
      double w = quad->weights[n] * area;
      for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
          // 直接使用排序好的 grad 数组，避免间接寻址
          double dot = grad[i][0] * grad[j][0] + grad[i][1] * grad[j][1];
          cc[i][j] += dot * w;
        }
      }
    }
  }

  /****************************************************************************
   * 棱边基函数与棱边基函数的点积
   ****************************************************************************/
  template <class ShapeFunc>
  void edgeBasDotedgeBas(const int face, const int order, const bool raw_order, double area,
                         const ShapeFunc *shapefunc, double *value,
                         const tbox::Array<tbox::Array<double> > &nodeGradBas) const {
    check(order);
    const int nbas = shapefunc->nbas();
    const int dim = shapefunc->dim();
    double(*bb)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    const Quad *quad = d_face_quad_table[order];
    std::fill(&(bb[0][0]), &(bb[0][0]) + nbas * nbas, 0.0);
    int faceorder = (raw_order) ? 1 : 0;
    for (int n = 0; n < quad->npoints; n++) {
      double phi[nbas][dim];
      shapefunc->basis(face, quad->points + n * quad->dim, &(phi[0][0]), nodeGradBas, faceorder);

      // 预乘面积权重
      double weight = quad->weights[n] * area;

      for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
          double v = 0.;
          for (int d = 0; d < dim; d++) {
            v += phi[i][d] * phi[j][d];
          }
          bb[i][j] += v * weight; // 直接累加
        }
      }
    }
  }

  /****************************************************************************
   * 棱边基函数旋度与棱边基函数旋度的点积
   ****************************************************************************/
  template <class ShapeFunc>
  void edgeCurlDotedgeCurl(const int face, const int order, const bool raw_order, double area,
                           const ShapeFunc *shapefunc, double *value,
                           const tbox::Array<tbox::Array<double> > &nodeGradBas) const {
    check(order);
    const int nbas = shapefunc->nbas();
    double(*cc)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    const Quad *quad = d_face_quad_table[order];
    std::fill(&(cc[0][0]), &(cc[0][0]) + nbas * nbas, 0.0);
    int faceorder = (raw_order) ? 1 : 0;
    for (int n = 0; n < quad->npoints; n++) {
      double curlphi[nbas][1];
      shapefunc->curl(face, quad->points + n * 3, &(curlphi[0][0]), nodeGradBas, faceorder);

      // 提前将面积与权重合并
      double weight = quad->weights[n] * area;

      for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
          double v = curlphi[i][0] * curlphi[j][0];
          cc[i][j] += v * weight;
        }
      }
    }
  }

  /****************************************************************************
   * 棱边基函数与结点基函数梯度的点积
   ****************************************************************************/
  template <class ShapeFunc>
  void edgeBasDotNodeGrad(const int face, const int order, const bool raw_order, double area,
                          const ShapeFunc *shapefunc, double *value,
                          const tbox::Array<tbox::Array<double> > &nodeGradBas) const {
    check(order);
    const int nbas = shapefunc->nbas();
    const int dim = shapefunc->dim();
    double(*bb)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    const Quad *quad = d_face_quad_table[order];
    std::fill(&(bb[0][0]), &(bb[0][0]) + nbas * nbas, 0.0);
    int faceorder = 0;
    int i_order[3] = {0, 2, 1};
    if (raw_order) {
      faceorder = 1;
      i_order[1] = 1;
      i_order[2] = 2;
    }
    // 预排序 nodeGradBas
    double grad[nbas][dim];
    for (int j = 0; j < nbas; j++) {
      grad[j][0] = nodeGradBas[i_order[j]][0];
      grad[j][1] = nodeGradBas[i_order[j]][1];
    }
    for (int n = 0; n < quad->npoints; n++) {
      double phi[nbas][dim];
      shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]), nodeGradBas, faceorder);
      double weight = quad->weights[n] * area;
      for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
          double v = 0.;
          for (int d = 0; d < dim; d++) {
            v += phi[i][d] * grad[j][d];
          }
          bb[i][j] += v * weight;
        }
      }
    }
  }

  /****************************************************************************
   * 结点基函数梯度与棱边基函数的点积
   ****************************************************************************/
  template <class ShapeFunc>
  void NodeGradDotedgeBas(const int face, const int order, const bool raw_order, double area,
                          const ShapeFunc *shapefunc, double *value,
                          const tbox::Array<tbox::Array<double> > &nodeGradBas) const {
    check(order);
    const int nbas = shapefunc->nbas();
    const int dim = shapefunc->dim();
    double(*bb)[nbas] = reinterpret_cast<double(*)[nbas]>(value);
    const Quad *quad = d_face_quad_table[order];
    std::fill(&(bb[0][0]), &(bb[0][0]) + nbas * nbas, 0.0);
    int faceorder = 0;
    int i_order[3] = {0, 2, 1};
    if (raw_order) {
      faceorder = 1;
      i_order[2] = 2;
      i_order[1] = 1;
    }
    double grad[nbas][dim];
    for (int i = 0; i < nbas; i++) {
      grad[i][0] = nodeGradBas[i_order[i]][0];
      grad[i][1] = nodeGradBas[i_order[i]][1];
    }
    for (int n = 0; n < quad->npoints; n++) {
      double phi[nbas][dim];
      shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]), nodeGradBas, faceorder);
      double weight = quad->weights[n] * area;

      for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
          double v = 0.;
          for (int d = 0; d < dim; d++) {
            v += grad[i][d] * phi[j][d];
          }
          bb[i][j] += v * weight;
        }
      }
    }
  }

  /****************************************************************************
   * 计算三角形面上的功率
   ****************************************************************************/
  template <class ShapeFunc, class TYPE, class TYPE2, class TYPE3>
  void faceEPower(const int face, const int order, const bool raw_order, double area,
                  const ShapeFunc *shapefunc, const TYPE *E_edge, const TYPE *E_node, TYPE2 *value,
                  const tbox::Array<tbox::Array<double> > &nodeGradBas, TYPE3 beta, double omega,
                  double mu) const {
    check(order);
    const int nbas = shapefunc->nbas();
    const int dim = shapefunc->dim();
    value[0] = TYPE2(0.0);
    const Quad *quad = d_face_quad_table[order];
    int faceorder = 0;
    int i_order[3] = {0, 2, 1};
    if (raw_order) {
      faceorder = 1;
      i_order[1] = 1;
      i_order[2] = 2;
    }
    double grad[nbas][dim];
    for (int i = 0; i < nbas; i++) {
      grad[i][0] = nodeGradBas[i_order[i]][0];
      grad[i][1] = nodeGradBas[i_order[i]][1];
    }
    const TYPE2 j_unit(0.0, 1.0);
    const TYPE2 j_omega_mu = j_unit * TYPE2(omega * mu);
    for (int n = 0; n < quad->npoints; n++) {
      double phi[nbas][dim]; // 现在 nbas 和 dim 是 const 常量，属于合法静态数组
      shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]), nodeGradBas, faceorder);
      // 横向场计算
      // et = beta * Et
      TYPE et[2] = {TYPE(0.0), TYPE(0.0)};
      for (int i = 0; i < nbas; i++) {
        et[0] += E_edge[i] * phi[i][0];
        et[1] += E_edge[i] * phi[i][1];
      }
      // 纵向场梯度计算
      // gradez = -j * Grad(Ez)
      TYPE gradez[2] = {TYPE(0.0), TYPE(0.0)};
      for (int i = 0; i < nbas; i++) {
        gradez[0] += E_node[i] * grad[i][0];
        gradez[1] += E_node[i] * grad[i][1];
      }

      TYPE2 Ht[2];
      Ht[0] = -(gradez[1] * j_unit + j_unit * et[1]) / j_omega_mu;
      Ht[1] = -(-gradez[0] * j_unit - j_unit * et[0]) / j_omega_mu;

      TYPE2 Et[2];
      Et[0] = et[0] / beta;
      Et[1] = et[1] / beta;
      // 电磁功率必须对磁场取复共轭
      TYPE2 S_flow = Et[0] * std::conj(Ht[1]) - Et[1] * std::conj(Ht[0]);
      value[0] += S_flow * TYPE2(quad->weights[n] * area);
    }
  }

  /////////////////////////////////////////模式功率（E叉乘H）/////////////////////////////////
  template <class ShapeFunc, class TYPE, class TYPE2, class TYPE3>
  void quadfaceEsquare(const int face, const int order, double area, const ShapeFunc *shapefunc,
                       const TYPE *E_edge, ////边上插值系数
                       const TYPE *E_node, /// 节点上插值系数
                       TYPE2 *value, tbox::Array<tbox::Array<double> > nodeGradBas, TYPE3 beta,
                       double omega, double mu) const {
    check(order);
    int nbas = shapefunc->nbas();
    int dim = shapefunc->dim();
    // TYPE(*bb) = (double(*)[1])value;
    TYPE2(*bb) = value;
    bzero(&(bb[0]), 1 * sizeof(*value));
    const Quad *quad = d_face_quad_table[order];
    /// 初始化bb
    bzero(&(bb[0]), 1 * sizeof(*value));
    int faceorder = 0; // 表示面的节点是否是逆时针编号
    int i_order[3] = {0, 2, 1};
    if (area > 0) {
      faceorder = 1;
      i_order[2] = 2;
      i_order[1] = 1;
    }

    for (int n = 0; n < quad->npoints; n++) {
      double phi[nbas][dim];
      shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]), nodeGradBas, faceorder);

      TYPE et[2] = {0, 0};          // 横向场//et=beta*Et
      for (int i = 0; i < 3; i++) { // 一个面上有三个基函数
        et[0] += E_edge[i] * phi[i][0];
        et[1] += E_edge[i] * phi[i][1];
      }
      //           cout<<"phi "<<phi[0][0]<<" "<<phi[0][1]<<endl;
      //           cout<<"phi "<<phi[1][0]<<" "<<phi[1][1]<<endl;
      //           cout<<"phi "<<phi[2][0]<<" "<<phi[2][1]<<endl;
      //           cout<<"Et "<<et[0]<<" "<<et[1]<<"   E_edge "<<E_edge[0]<<"
      //           "<<E_edge[1]<<" "<<E_edge[2]<<endl;
      TYPE gradez[2] = {0, 0}; // 纵向场ez=-j*Ez
      for (int i = 0; i < 3; i++) {
        gradez[0] +=
            E_node[i] *
            nodeGradBas
                [i_order[i]]
                [0]; /// 节点2，3在area为负值时交换了顺序，故求出梯度顺序是0，2，1不同于面单元所规定的顺序0，1，2
        gradez[1] += E_node[i] * nodeGradBas[i_order[i]][1];
      }
      //           cout<<E_node[0]<<" "<<E_node[1]<<E_node[2]<<endl;
      TYPE2 v = 0.;
      TYPE2 Ht[2] = {0,
                     0}; // Hx=(dEz/dy+j*beta*Ety)/(j*omega*mu),Hy=(-dEz/dx-j*beta*Etx)/(j*omega*mu)
      Ht[0] =
          -(gradez[1] * dcomplex(0, 1) + dcomplex(0, 1) * et[1]) / (dcomplex(0, 1) * omega * mu);
      Ht[1] =
          -(-gradez[0] * dcomplex(0, 1) - dcomplex(0, 1) * et[0]) / (dcomplex(0, 1) * omega * mu);
      //            cout<<Ht[0]<<" "<<Ht[1]<<endl;
      TYPE2 Et[2] = {0, 0};
      Et[0] = et[0] / beta;
      Et[1] = et[1] / beta;
      // cout<<"Et "<<Et[0]<<" "<<Et[1]<<endl;
      TYPE2 S_flow = 0; // E叉乘H，入射功率
      S_flow = Et[0] * Ht[1] - Et[1] * Ht[0];
      // cout<<"yita"<<Et[0]/Ht[1]<<"  "<<-Et[1]/Ht[0]<<endl;
      v += S_flow;
      // cout<<"Enorm "<<sqrt(v)<<endl;
      v *= quad->weights[n];
      (*bb) = (*bb) + v;
    }
    (*bb) = (*bb) * abs(area);
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
    if (order > max_order)
      TBOX_ERROR("Quadrature rules of order " << order << " for tetrahedra not implememted.\n");
  }

  const hier::Patch<NDIM> &d_patch;
  tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_volume;
  tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;

  tbox::Array<int> can_ext, can_idx;
  tbox::Array<int> fan_ext, fan_idx;
  tbox::Array<int> caf_ext, caf_idx;
  tbox::Array<int> ean_ext, ean_idx;
  tbox::Array<int> cae_ext, cae_idx;

  Quad **d_face_quad_table;
  Quad **d_edge_quad_table;
};
} // namespace appu
} // namespace JAUMIN

#endif // included_appu_TetQuad
