#pragma once

#include <assert.h>
#include <math.h>

/**
 * @brief 求向量vec1和vec2的内积。
 * @param n
 * @param vec1
 * @param vec2
 * @return
 */
template <class TYPE>
static TYPE dotProduct(int n,
        const TYPE* vec1, const TYPE* vec2)
{
    TYPE sum = TYPE(0);
    for (int i = 0; i < n; i++)
        sum += vec1[i] * vec2[i];
    return sum;
}
template <class TYPE, class TYPE2>
static TYPE dotProduct(int n,
        const TYPE* vec1, const TYPE2* vec2)
{
    TYPE sum = TYPE(0);
    for (int i = 0; i < n; i++)
        sum += vec1[i] * TYPE(vec2[i]);
    return sum;
}

/**
 * @brief 求向量vec1和vec2的叉积。
 * @param n
 * @param vec1
 * @param vec2
 * @param vec
 */
template <class TYPE>
static void crossProduct(int n,
        const TYPE* vec1, const TYPE* vec2,
        TYPE* vec)
{
    assert(3 == n);
    for (int k = 0; k < n; k++)
    {
        int i = (k + 1) % n;
        int j = (k + 2) % n;
        vec[k] = vec1[i] * vec2[j] - vec1[j] * vec2[i];
    }
}
template <class TYPE, class TYPE2>
static void crossProduct(int n,
        const TYPE* vec1, const TYPE2* vec2,
        TYPE* vec)
{
    assert(3 == n);
    for (int k = 0; k < n; k++)
    {
        int i = (k + 1) % n;
        int j = (k + 2) % n;
        vec[k] = vec1[i] * TYPE(vec2[j]) - vec1[j] * TYPE(vec2[i]);
    }
}
/**
 * @brief 取出向量的方向
 * @param n 向量的长度
 * @param data 向量的内容
 * @return 返回true，如果x>0 || (x==0 && y>0) || (x==0 && y==0 && z>0)；
 * 否则返回false
 */
template<class TYPE>
bool vectorDirection(int n, const TYPE* data) {
    for (int i = 0; i < n; i++) {
        if (data[i] > 0)
            return true;
        else if (data[i] == 0)
            continue;
        else
            return false;
    }
    assert(false);
    return false;
}

/**
 * @brief 将向量的内容反向存储
 * @param n 向量长度
 * @param data 向量内容
 */
template<class TYPE>
void reverseOrder(int n, TYPE* data) {
    if (n <= 1) return;
    for (int i = 0; i < n / 2; i++) {
        int j = n - 1 - i;
        assert(i < j);
        TYPE tmp = data[i];
        data[i] = data[j];
        data[j] = tmp;
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
void outerNormal(
        int ncellnode, const int* cellnodes,
        int nfacenode, const int* facenodes,
        int nnodes, const double* nodecoords,
        int ndim, double* normal, double* normal_raw = 0);


/**
 * @brief barycentricToSpace 将barycentric中的重心坐标（ndim+1维）转为空间坐标，
 * 空间中的4个点由nodes表示。
 * @param ndim
 * @param nnodes
 * @param nodecoords
 * @param nodes
 * @param barycentric
 * @param spacecoords
 */
void barycentricToSpace(int ndim, int nnodes, const double* nodecoords,
        const int* spacenodes, const double* barycentric, double* spacecoords);
