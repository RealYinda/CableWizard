#include "JAUMIN_Macros.h"
#include "TriGeom.h"

/*
在端口面上增加一个二维坐标系，其中的u与正交坐标系对齐，nuv构成右手坐标系
*/
void Transfer3DCoordTo2D(const tbox::Array<hier::DoubleVector<NDIM> >& coord_3d,
                         const tbox::Array<double>& direction,
                         tbox::Array<hier::DoubleVector<NDIM - 1> >& coord_2d) {
  // 获取并归一化端口的外法向向量
  double n[3] = {direction[0], direction[1], direction[2]};
  NORMALIZATION(n);
  double u[3] = {0}, v[3] = {0};
  
  // 选取参考平面
  if (fabs(n[0]) < 0.9) {
    u[0] = 1; u[1] = 0; u[2] = 0; // 参考 X 轴
  } else {
    u[0] = 0; u[1] = 1; u[2] = 0; // 参考 Y 轴 
  }
  
  // Gram-Schmidt 正交化：去掉 u 在 n 上的投影，使 u 垂直于 n
  double dot_un = u[0] * n[0] + u[1] * n[1] + u[2] * n[2];
  u[0] -= dot_un * n[0];
  u[1] -= dot_un * n[1];
  u[2] -= dot_un * n[2];
  NORMALIZATION(u);
  
  // 利用叉乘生成 v = n x u
  CROSS_PRODUCT(n, u, v);

  int nnode = 3;
  for (int i = 0; i < nnode; i++) {
    coord_2d[i][0] = coord_3d[i][0] * u[0] + coord_3d[i][1] * u[1] + coord_3d[i][2] * u[2];
    coord_2d[i][1] = coord_3d[i][0] * v[0] + coord_3d[i][1] * v[1] + coord_3d[i][2] * v[2];
  }
}

void GradientOn2DCoord(const tbox::Array<hier::DoubleVector<NDIM - 1> >& coord_2d,
                       tbox::Array<tbox::Array<double> >& gradient) {
  tbox::Array<hier::DoubleVector<NDIM - 1> > real_vertex = coord_2d;

  double area = AREA2DQUAD(real_vertex[0], real_vertex[1], real_vertex[2]);

  tbox::Array<double> val_0(2);
  GRADIENT2D(val_0, real_vertex[1], real_vertex[2], area);
  gradient[0] = val_0;

  tbox::Array<double> val_1(2);
  GRADIENT2D(val_1, real_vertex[2], real_vertex[0], area);
  gradient[1] = val_1;

  tbox::Array<double> val_2(2);
  GRADIENT2D(val_2, real_vertex[0], real_vertex[1], area);
  gradient[2] = val_2;
}