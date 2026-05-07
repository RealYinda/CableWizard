#include "GridInfo.h"

/**
 * @brief 求面的右手法向，索引都是从0开始
 * @param nfacenode 面中点的个数
 * @param facenodes 面中的点的在网格中的编号
 * @param nnodes 网格节点个数
 * @param nodecoords 网格节点坐标
 * @param dim 问题维数
 * @param normal [out] 输出面的右手法向
 */
void rightHandNormal(
        int nfacenode, const int* facenodes,
        int nnodes, const double* nodecoords,
        int ndim, double* normal, double* normal_raw)
{
    assert(ndim == 3);
    double* vec1 = new double[ndim];
    double* vec2 = new double[ndim];
    for (int i = 0; i<ndim; i++)
    {
        if (i < ndim)
        {
            vec1[i] = nodecoords[facenodes[1]*ndim+i]
                    - nodecoords[facenodes[0]*ndim+i];
            vec2[i] = nodecoords[facenodes[2]*ndim+i]
                    - nodecoords[facenodes[1]*ndim+i];
        }
        else
        {
            vec1[i] = 0.0;
            vec2[i] = 0.0;
        }
    }
    crossProduct(ndim, vec1, vec2, normal);
    // 归一化，将长度化为1
    double sum = dotProduct(ndim, normal, normal);
    for (int i=0; i<ndim; i++) {
        if (normal_raw) normal_raw[i] = normal[i];
        normal[i] /= sqrt(sum);
    }
    delete[] vec1;
    delete[] vec2;
}

/**
 * @brief 取出单元中的一条边，一个点位于给定的面上，另一个点不在面上。
 * @param ncellnode 单元中点的个数
 * @param cellnodes 单元中的点在网格中的编号
 * @param nfacenode 面中点的个数
 * @param facenodes 面中的点的在网格中的编号
 * @param inface [out] 输出位于面上点的编号
 * @param outface [out] 输出位于面外的点的编号
 */
void getLineInCell(
        int ncellnode, const int* cellnodes,
        int nfacenode, const int* facenodes,
        int* inface, int* outface)
{
    int* facenode_indicator = new int[ncellnode];
    for (int j=0; j<ncellnode; j++)
        facenode_indicator[j] = 0;
    /*
     * 将单元中属于面的点标记出来。
     */
    for (int i=0; i<nfacenode; i++)
    {
        bool is_in_cell = false;
        for (int j=0; j<ncellnode; j++)
        {
            if (cellnodes[j] == facenodes[i])
            {
                facenode_indicator[j] = 1;
                is_in_cell = true;
                break;
            }
        }
        assert(is_in_cell);
    }
#ifdef TEST_IN_IMPLIMENTATION
    int cnt = 0;
    for (j=0; j<ncellnode; j++)
        cnt += facenode_indicator[j];
    assert(cnt == nfacenode);
#endif
    int j;
    /*
     * 取出一个属于面上的点。
     */
    for (j=0; j<ncellnode; j++)
    {
        if (facenode_indicator[j] == 1)
        {
            *inface = cellnodes[j];
            break;
        }
    }
    assert(j < ncellnode);
    /*
     * 取出一个不属于面上的点。
     */
    for (j=0; j<ncellnode; j++)
    {
        if (facenode_indicator[j] == 0)
        {
            *outface = cellnodes[j];
            break;
        }
    }
    assert(j < ncellnode);
    delete[] facenode_indicator;
}

/**
 * @brief 取单元中的一条边，索引都是从0开始
 * @param ncellnode 单元中点的个数
 * @param cellnodes 单元中的点在网格中的编号
 * @param nfacenode 面中点的个数
 * @param facenodes 面中的点的在网格中的编号
 * @param nnodes 网格节点个数
 * @param nodecoords 网格节点坐标
 * @param dim 问题维数
 * @param dir [out] 单元上边的方向, 与这个面的外法向方向相同
 */
void getLineInCell(
        int ncellnode, const int* cellnodes,
        int nfacenode, const int* facenodes,
        int nnodes, const double* nodecoords,
        int ndim, double* dir)
{
    // 计算出单元上的一条边，其中一个点位于面上，一个点位于面外。
    int inface, outface;
    getLineInCell(ncellnode, cellnodes,
            nfacenode, facenodes,
            &inface, &outface);
    for (int i=0; i<ndim; i++)
    {
        dir[i] = nodecoords[inface*ndim+i]
               - nodecoords[outface*ndim+i];
    }
}

/**
 * @brief 求单元中面的外法向，索引都是从0开始
 * @param ncellnode 单元中点的个数
 * @param cellnodes 单元中的点在网格中的编号
 * @param nfacenode 面中点的个数
 * @param facenodes 面中的点的在网格中的编号
 * @param nnodes 网格节点个数
 * @param nodecoords 网格节点坐标
 * @param dim 问题维数
 * @param normal [out] 输出面的外法向
 */
void outerNormal(int ncellnode, const int* cellnodes,
        int nfacenode, const int* facenodes,
        int nnodes, const double* nodecoords,
        int ndim, double* normal, double *normal_raw)
{
    rightHandNormal(nfacenode, facenodes, nnodes, nodecoords, ndim, normal, normal_raw);
    // 计算出单元上的一条边，其中一个点位于面上，一个点位于面外。
    double* dir = new double[ndim];
    getLineInCell(ncellnode, cellnodes, nfacenode, facenodes,
            nnodes, nodecoords, ndim, dir);
    // 如果右手法向与边的方向小于90度，则右手法向是外法向。
    if (dotProduct(ndim, normal, dir) < 0)
    {
        for (int i=0; i<ndim; i++)
        {
            normal[i] = -normal[i];
        }
    }
    delete[] dir;
}

/**
 * @brief 求单元中面的全局法向，从外法向得到的。
 * 如果外法向的方向为正，则不改变，否则取反方向。
 * 索引都是从0开始。
 * @param ncellnode 单元中点的个数
 * @param cellnodes 单元中的点在网格中的编号
 * @param nfacenode 面中点的个数
 * @param facenodes 面中的点的在网格中的编号
 * @param nnodes 网格节点个数
 * @param nodecoords 网格节点坐标
 * @param dim 问题维数
 * @param normal [out] 输出面的外法向
 */
void globalNormal(int ncellnode, const int* cellnodes, int nfacenode,
        const int* facenodes, int nnodes, const double* nodecoords, int ndim,
        double* normal) {
    outerNormal(ncellnode, cellnodes, nfacenode, facenodes, nnodes, nodecoords,
            ndim, normal);
    if (!vectorDirection(ndim, normal)) {
        for (int i = 0; i < ndim; i++) {
            normal[i] = -normal[i];
        }
    }
}


void barycentricToSpace(
        int ndim,
        int nnodes,
        const double *nodecoords,
        const int *spacenodes,
        const double *barycentric,
        double *spacecoords)
{
    for (int i=0; i<ndim; i++) {
        spacecoords[i] = 0;
    }

    for (int i=0; i<ndim+1; i++) {
        assert(spacenodes[i] < nnodes);
        const double* space_i = nodecoords + ndim * spacenodes[i];
        for (int j=0; j<ndim; j++) {
            spacecoords[j] += barycentric[i] * space_i[j];
        }
    }
}
