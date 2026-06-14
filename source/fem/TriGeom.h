#ifndef TRIGEOM_H
#define TRIGEOM_H
#include "TriQuad.h"
using namespace JAUMIN;

#define AREA2DQUAD(a, b, c)                                                    \
  ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))

#define GRADIENT2D(val, a, b, area2D)                                          \
  do {                                                                         \
    (val)[0] = ((a)[1] - (b)[1]) / (area2D);                                   \
    (val)[1] = ((b)[0] - (a)[0]) / (area2D);                                   \
  } while (0)

void Transfer3DCoordTo2D(const tbox::Array<hier::DoubleVector<NDIM> >& coord_3d,
                         const tbox::Array<double>& direction,
                         tbox::Array<hier::DoubleVector<NDIM - 1> >& coord_2d);

void GradientOn2DCoord(const tbox::Array<hier::DoubleVector<NDIM - 1> >& coord_2d,
                       tbox::Array<tbox::Array<double> >& gradient);
#endif