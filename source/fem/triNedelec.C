//
// 文件名:      Nedelec.C
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

#include "triNedelec.h"

namespace JAUMIN {
namespace appu {

triNedelec::triNedelec(const hier::Patch<NDIM>& patch,
                 tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
                 tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian)
    : d_patch(patch), d_edge_order(edge_order), d_cell_jacobian(cell_jacobian) {
  patch.getPatchTopology()->getFaceAdjacencyEdges(d_face_edge_ext,
                                                  d_face_edge_idx);
  patch.getPatchTopology()->getEdgeAdjacencyNodes(d_edge_node_ext,
                                                  d_edge_node_idx);
  patch.getPatchTopology()->getFaceAdjacencyNodes(d_face_node_ext,
                                                  d_face_node_idx);
  /// 取出本地Patch的结点坐标数组.
  d_node_coord =patch.getPatchGeometry()->getNodeCoordinates();

}

int triNedelec::nbas() const { return 3; }

int triNedelec::dim() const { return 2; }

void triNedelec::basis(const int face, const double* lambda, double* values,const tbox::Array<tbox::Array<double> > &nodeGradBas,int faceorder) const
/*
 * \phi_ij =  \lambda_i \grad \lambda_j - \lambda_j \grad \lambda_i
 *
 * zgd： lambda_i,高斯积分点处基函数i的值。
 */
{
//  double(*nabla)[NDIM + 1] =
//      (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));
  double(*phi)[2] = (double(*)[2])values;

  int nedge = 3;
  for (int ie = 0; ie < nedge; ie++) {
    int edge = d_face_edge_idx[d_face_edge_ext[face] + ie];//单元的第几条边
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];//边的一个节点的全局编号
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];//边的一个节点的全局编号
//    cout<<v0<<" "<<v1<<endl;
//    double x0, y0, z0, x1, y1, z1;
//    x0 = (*d_node_coord)(0, v0);
//    y0 = (*d_node_coord)(1, v0);
//    z0 = (*d_node_coord)(2, v0);
//    x1 = (*d_node_coord)(0, v1);
//    y1 = (*d_node_coord)(1, v1);
//    z1 = (*d_node_coord)(2, v1);
    //double l=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
    bool in_order = (*d_edge_order)(0, edge);
    if (!in_order) {
      int t = v0;
      v0 = v1;
      v1 = t;//交换顺序
    }
    //cout<<v0<<" "<<v1<<endl;
    int lv0, lv1;
    lv0 = lv1 = -1;
    //确定边的一个节点的局部编号
    for (int n = 0; n < 3; n++) {
      int nidx = d_face_node_idx[d_face_node_ext[face] + n];//单元的节点
      int n_order[3]={0,1,2};
      if(!faceorder)
      {
         n_order[2]=1;
         n_order[1]=2;
      }
     //cout<<nidx<<" ";
      if (nidx == v0)
            lv0 = n_order[n];
      if (nidx == v1)
            lv1 = n_order[n];
    }
    //cout<<endl;
    phi[ie][0] = (lambda[lv0] * nodeGradBas[lv1][0] - lambda[lv1] * nodeGradBas[lv0][0]);//lanmda[]表示节点基函数在高斯积分点的值
    phi[ie][1] = (lambda[lv0] * nodeGradBas[lv1][1] - lambda[lv1] * nodeGradBas[lv0][1]);//nabla表示基函数的梯度在高斯积分点的值
//    cout<<"ie "<<ie<<" "<<v0<<" "<<v1<<" "<<lv0<<" "<<lv1<<"|| "
//                   <<nodeGradBas[lv0][0]<<" "<<nodeGradBas[lv0][1]<<" "
//                  <<nodeGradBas[lv1][0]<<" "<<nodeGradBas[lv1][1]<<endl;
  }
//  cout<<nodeGradBas[0][0]<<" "<<nodeGradBas[0][1]<<endl;
//  cout<<nodeGradBas[1][0]<<" "<<nodeGradBas[1][1]<<endl;
//   cout<<nodeGradBas[2][0]<<" "<<nodeGradBas[2][1]<<endl;
}

void triNedelec::curl(const int face, const double* lambda, double* values,const tbox::Array<tbox::Array<double> > &nodeGradBas,int faceorder) const
/*
 * \curl \phi_ij =  2\grad\lambda_i x \grad\lambda_j
 */
{
  NULL_USE(lambda);
//  double(*nabla)[NDIM + 1] =
//      (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, face)));
  double(*curlphi)[1] = (double(*)[1])values;
  int nedge = 3;
  for (int ie = 0; ie < nedge; ie++) {
    int edge = d_face_edge_idx[d_face_edge_ext[face] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];
//    double x0, y0, z0, x1, y1, z1;
//    x0 = (*d_node_coord)(0, v0);
//    y0 = (*d_node_coord)(1, v0);
//    z0 = (*d_node_coord)(2, v0);
//    x1 = (*d_node_coord)(0, v1);
//    y1 = (*d_node_coord)(1, v1);
//    z1 = (*d_node_coord)(2, v1);
    //double l=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
    bool in_order = (*d_edge_order)(0, edge);
    if (!in_order) {
      int t = v0;
      v0 = v1;
      v1 = t;
    }

    int lv0, lv1;
    lv0 = lv1 = -1;
    for (int n = 0; n < 3; n++) {
      int nidx = d_face_node_idx[d_face_node_ext[face] + n];
      int n_order[3]={0,1,2};
      if(!faceorder)
      {
         n_order[2]=1;
         n_order[1]=2;
      }
      //cout<<nidx<<" ";
       if (nidx == v0) lv0 = n_order[n];
       if (nidx == v1) lv1 = n_order[n];
    }

    double a = nodeGradBas[lv0][0];//1,1
    double b = nodeGradBas[lv0][1];//1,2
    double d = nodeGradBas[lv1][0];//2,1
    double e = nodeGradBas[lv1][1];//2,2

    curlphi[ie][0] = (a * e - d * b) *2;//1,1*2,2-2,1*1,2

  //cout<<"ie "<<ie<<" "<<v0<<" "<<v1<<" "<<lv0<<" "<<lv1<<"|| "<<curlphi[ie][0]<<endl;
  }

}

}  // namespace appu
}  // namespace JAUMIN
