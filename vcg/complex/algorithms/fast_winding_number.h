/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef __VCG_TRI_FAST_WINDING_NUMBER
#define __VCG_TRI_FAST_WINDING_NUMBER

#include <vcg/space/index/aabb_binary_tree/aabb_binary_tree.h>
#include <cmath>
#include <vector>
#include <algorithm>

#ifndef FAST_WINDING_OMP_MIN_VALUE
#define FAST_WINDING_OMP_MIN_VALUE 1000
#endif

#ifndef TAYLOR_SERIES_ORDER
#define TAYLOR_SERIES_ORDER 2
#endif

namespace vcg {
namespace tri {


/**
 * BoxData stores precomputed values for each node in the hierarchy
 * Following the structure from UT_SolidAngle.cpp
 */
template<typename T, typename S>
struct BoxData
{
    void clear()
    {
        memset(this, 0, sizeof(*this));
    }

    // An upper bound on the squared distance from averageP to the farthest point in the box
    S maxPDist2;

    // Centre of mass of the mesh surface in this box
    Point3<T> averageP;

    // Unnormalized, area-weighted normal of the mesh in this box
    Point3<T> N;

#if TAYLOR_SERIES_ORDER >= 1
    // Values for Omega_1 (first-order Taylor expansion)
    Point3<T> NijDiag;  // Nxx, Nyy, Nzz
    T Nxy_Nyx;          // Nxy+Nyx
    T Nyz_Nzy;          // Nyz+Nzy
    T Nzx_Nxz;          // Nzx+Nxz
#endif

#if TAYLOR_SERIES_ORDER >= 2
    // Values for Omega_2 (second-order Taylor expansion)
    Point3<T> NijkDiag;     // Nxxx, Nyyy, Nzzz
    T sumPermuteNxyz;       // (Nxyz+Nxzy+Nyzx+Nyxz+Nzxy+Nzyx) = 2*(Nxyz+Nyzx+Nzxy)
    T N2xxy_Nyxx;           // 2*Nxxy+Nyxx
    T N2xxz_Nzxx;           // 2*Nxxz+Nzxx
    T N2yyz_Nzyy;           // 2*Nyyz+Nzyy
    T N2yyx_Nxyy;           // 2*Nyyx+Nxyy
    T N2zzx_Nxzz;           // 2*Nzzx+Nxzz
    T N2zzy_Nyzz;           // 2*Nzzy+Nyzz
#endif

    // Additional data needed during tree construction
    T area;                 // Total area in this box
    Point3<T> areaP;        // Area-weighted position
};


/**
 * Fast Winding Number class
 * Implements hierarchical evaluation of the generalized winding number
 * using a BVH tree with Taylor series approximations
 */
template<typename MeshType>
class FastWindingNumber
{
public:
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef BoxData<ScalarType, ScalarType> BoxDataType;
    typedef AABBBinaryTree<FaceType, ScalarType, BoxDataType> TreeType;
    typedef typename TreeType::NodeType NodeType;

    FastWindingNumber();
    ~FastWindingNumber();

    // Initialize the hierarchical structure
    void init(const MeshType &mesh, const int order = 2);

    // Compute winding number for a query point
    ScalarType computeSolidAngle(const CoordType &queryPoint, const ScalarType accuracyScale = 2.0) const;

    // Clear the structure
    void clear();

private:
    TreeType tree;
    const MeshType *pMesh;
    int myOrder;
    int nTriangles;

    // Precompute coefficients for a leaf node (single triangle)
    void precomputeLeaf(NodeType *node, const FaceType *face);

    // Precompute coefficients for an internal node
    void precomputeInternal(NodeType *node);

    // Recursively precompute all nodes
    void precomputeNode(NodeType *node);

    // Compute solid angle for a single triangle (exact calculation)
    ScalarType computeTriangleSolidAngle(const CoordType &queryPoint, const FaceType *face) const;

    // Evaluate using Taylor series approximation
    ScalarType evaluateApproximation(const CoordType &queryPoint, const NodeType *node) const;

    // Traverse the tree and accumulate solid angles
    ScalarType traverseTree(const CoordType &queryPoint, const NodeType *node, const ScalarType accuracyScale2) const;
};


// ============================================================================
// Implementation
// ============================================================================

template<typename MeshType>
FastWindingNumber<MeshType>::FastWindingNumber()
    : pMesh(nullptr)
    , myOrder(2)
    , nTriangles(0)
{
}

template<typename MeshType>
FastWindingNumber<MeshType>::~FastWindingNumber()
{
}

template<typename MeshType>
void FastWindingNumber<MeshType>::clear()
{
    tree.Clear();
    pMesh = nullptr;
    myOrder = 2;
    nTriangles = 0;
}

// Compute solid angle of a single triangle using Van Oosterom & Strackee formula
template<typename MeshType>
typename FastWindingNumber<MeshType>::ScalarType 
FastWindingNumber<MeshType>::computeTriangleSolidAngle(const CoordType &queryPoint, const FaceType *face) const
{
    if (face->IsD()) return 0.0;

    CoordType v0 = face->cP(0);
    CoordType v1 = face->cP(1);
    CoordType v2 = face->cP(2);

    CoordType a = v0 - queryPoint;
    CoordType b = v1 - queryPoint;
    CoordType c = v2 - queryPoint;

    ScalarType la = a.Norm();
    ScalarType lb = b.Norm();
    ScalarType lc = c.Norm();

    if (la < 1e-10 || lb < 1e-10 || lc < 1e-10)
        return 0.0;

    ScalarType det = a * (b ^ c);
    ScalarType dp_bc = b * c;
    ScalarType dp_ca = c * a;
    ScalarType dp_ab = a * b;

    ScalarType denominator = la * lb * lc + la * dp_bc + lb * dp_ca + lc * dp_ab;
    ScalarType solidAngle = 2.0 * atan2(det, denominator);

    return solidAngle;
}

// Precompute data for a leaf node (single triangle)
template<typename MeshType>
void FastWindingNumber<MeshType>::precomputeLeaf(NodeType *node, const FaceType *face)
{
    if (!node || !face || face->IsD())
        return;

    BoxDataType &data = node->auxData;
    data.clear();

    // Get triangle vertices
    CoordType v0 = face->cP(0);
    CoordType v1 = face->cP(1);
    CoordType v2 = face->cP(2);

    // Compute triangle area and normal
    CoordType edge01 = v1 - v0;
    CoordType edge02 = v2 - v0;
    CoordType crossProd = edge01 ^ edge02;
    ScalarType area2 = crossProd.Norm();
    ScalarType area = area2 * 0.5;

    if (area < 1e-10)
    {
        data.area = 0.0;
        return;
    }

    // Area-weighted normal (unnormalized)
    data.N = crossProd * 0.5;
    data.area = area;

    // Barycenter
    CoordType barycenter = (v0 + v1 + v2) / 3.0;
    data.averageP = barycenter;
    data.areaP = barycenter * area;

    // Compute max distance squared
    ScalarType d0 = (v0 - barycenter).SquaredNorm();
    ScalarType d1 = (v1 - barycenter).SquaredNorm();
    ScalarType d2 = (v2 - barycenter).SquaredNorm();
    data.maxPDist2 = std::max(d0, std::max(d1, d2));

#if TAYLOR_SERIES_ORDER >= 1
    // Compute second-order moments for Omega_1
    CoordType p = barycenter;
    ScalarType nx = data.N[0], ny = data.N[1], nz = data.N[2];
    
    data.NijDiag[0] = nx * p[0];  // Nxx
    data.NijDiag[1] = ny * p[1];  // Nyy
    data.NijDiag[2] = nz * p[2];  // Nzz
    
    data.Nxy_Nyx = nx * p[1] + ny * p[0];
    data.Nyz_Nzy = ny * p[2] + nz * p[1];
    data.Nzx_Nxz = nz * p[0] + nx * p[2];
#endif

#if TAYLOR_SERIES_ORDER >= 2
    // Compute third-order moments for Omega_2
    ScalarType px = p[0], py = p[1], pz = p[2];
    
    data.NijkDiag[0] = nx * px * px;  // Nxxx
    data.NijkDiag[1] = ny * py * py;  // Nyyy
    data.NijkDiag[2] = nz * pz * pz;  // Nzzz
    
    data.sumPermuteNxyz = 2.0 * (nx * py * pz + ny * pz * px + nz * px * py);
    data.N2xxy_Nyxx = 2.0 * nx * px * py + ny * px * px;
    data.N2xxz_Nzxx = 2.0 * nx * px * pz + nz * px * px;
    data.N2yyz_Nzyy = 2.0 * ny * py * pz + nz * py * py;
    data.N2yyx_Nxyy = 2.0 * ny * py * px + nx * py * py;
    data.N2zzx_Nxzz = 2.0 * nz * pz * px + nx * pz * pz;
    data.N2zzy_Nyzz = 2.0 * nz * pz * py + ny * pz * pz;
#endif
}

// Precompute data for an internal node by combining children
template<typename MeshType>
void FastWindingNumber<MeshType>::precomputeInternal(NodeType *node)
{
    if (!node || node->IsLeaf())
        return;

    BoxDataType &data = node->auxData;
    data.clear();

    // First, sum zero-order moments to find parent center
    for (int i = 0; i < 2; ++i)
    {
        NodeType *child = node->children[i];
        if (!child) continue;

        BoxDataType &childData = child->auxData;

        data.N += childData.N;
        data.area += childData.area;
        data.areaP += childData.areaP;
    }

    // Compute average position for parent
    if (data.area > 1e-10)
    {
        data.averageP = data.areaP / data.area;
    }

#if TAYLOR_SERIES_ORDER >= 1
    // Now adjust higher-order moments for displacement from child centers to parent center
    for (int i = 0; i < 2; ++i)
    {
        NodeType *child = node->children[i];
        if (!child) continue;

        BoxDataType &childData = child->auxData;
        
        // Displacement from child center to parent center
        CoordType displacement = childData.averageP - data.averageP;
        CoordType N = childData.N;
        ScalarType dx = displacement[0], dy = displacement[1], dz = displacement[2];
        ScalarType Nx = N[0], Ny = N[1], Nz = N[2];
        
        // Add child Nij and adjust for displacement
        // Nij_diag(parent) = Nij_diag(child) + Ni * dj (for diagonal terms i=j)
        data.NijDiag[0] += childData.NijDiag[0] + Nx * dx;
        data.NijDiag[1] += childData.NijDiag[1] + Ny * dy;
        data.NijDiag[2] += childData.NijDiag[2] + Nz * dz;
        
        // For off-diagonal terms
        ScalarType Nxy = Nx * dy;
        ScalarType Nyx = Ny * dx;
        ScalarType Nyz = Ny * dz;
        ScalarType Nzy = Nz * dy;
        ScalarType Nzx = Nz * dx;
        ScalarType Nxz = Nx * dz;
        
        data.Nxy_Nyx += childData.Nxy_Nyx + Nxy + Nyx;
        data.Nyz_Nzy += childData.Nyz_Nzy + Nyz + Nzy;
        data.Nzx_Nxz += childData.Nzx_Nxz + Nzx + Nxz;

#if TAYLOR_SERIES_ORDER >= 2
        // Adjust Nijk for the change in center P
        // Nijk_diag = Nijk_child + 2*di*Nij_child + di*di*Ni
        data.NijkDiag[0] += childData.NijkDiag[0] + 2.0 * dx * childData.NijDiag[0] + dx * dx * Nx;
        data.NijkDiag[1] += childData.NijkDiag[1] + 2.0 * dy * childData.NijDiag[1] + dy * dy * Ny;
        data.NijkDiag[2] += childData.NijkDiag[2] + 2.0 * dz * childData.NijDiag[2] + dz * dz * Nz;
        
        data.sumPermuteNxyz += childData.sumPermuteNxyz +
            dx * (Nyz + Nzy + childData.Nyz_Nzy) +
            dy * (Nzx + Nxz + childData.Nzx_Nxz) +
            dz * (Nxy + Nyx + childData.Nxy_Nyx);
        
        data.N2xxy_Nyxx += childData.N2xxy_Nyxx +
            2.0 * (dy * childData.NijDiag[0] + dx * Nxy + Nx * dx * dy) +
            2.0 * Nyx * dx + Ny * dx * dx;
            
        data.N2xxz_Nzxx += childData.N2xxz_Nzxx +
            2.0 * (dz * childData.NijDiag[0] + dx * Nxz + Nx * dx * dz) +
            2.0 * Nzx * dx + Nz * dx * dx;
            
        data.N2yyz_Nzyy += childData.N2yyz_Nzyy +
            2.0 * (dz * childData.NijDiag[1] + dy * Nyz + Ny * dy * dz) +
            2.0 * Nzy * dy + Nz * dy * dy;
            
        data.N2yyx_Nxyy += childData.N2yyx_Nxyy +
            2.0 * (dx * childData.NijDiag[1] + dy * Nyx + Ny * dy * dx) +
            2.0 * Nxy * dy + Nx * dy * dy;
            
        data.N2zzx_Nxzz += childData.N2zzx_Nxzz +
            2.0 * (dx * childData.NijDiag[2] + dz * Nzx + Nz * dz * dx) +
            2.0 * Nxz * dz + Nx * dz * dz;
            
        data.N2zzy_Nyzz += childData.N2zzy_Nyzz +
            2.0 * (dy * childData.NijDiag[2] + dz * Nzy + Nz * dz * dy) +
            2.0 * Nyz * dz + Ny * dz * dz;
#endif
    }
#endif

    // Compute max distance for internal node
    CoordType boxMin(node->boxCenter[0] - node->boxHalfDims[0],
                     node->boxCenter[1] - node->boxHalfDims[1],
                     node->boxCenter[2] - node->boxHalfDims[2]);
    CoordType boxMax(node->boxCenter[0] + node->boxHalfDims[0],
                     node->boxCenter[1] + node->boxHalfDims[1],
                     node->boxCenter[2] + node->boxHalfDims[2]);

    // Check all 8 corners of the box
    ScalarType maxDist2 = 0.0;
    for (int ix = 0; ix < 2; ++ix)
        for (int iy = 0; iy < 2; ++iy)
            for (int iz = 0; iz < 2; ++iz)
            {
                CoordType corner(ix ? boxMax[0] : boxMin[0],
                                 iy ? boxMax[1] : boxMin[1],
                                 iz ? boxMax[2] : boxMin[2]);
                ScalarType dist2 = (corner - data.averageP).SquaredNorm();
                maxDist2 = std::max(maxDist2, dist2);
            }

    data.maxPDist2 = maxDist2;
}


// Recursively precompute all nodes (post-order traversal)
template<typename MeshType>
void FastWindingNumber<MeshType>::precomputeNode(NodeType *node)
{
    if (!node)
        return;

    if (node->IsLeaf())
    {
        // Leaf node - process and accumulate all triangles
        BoxDataType &data = node->auxData;
        data.clear();
        
        std::vector<CoordType> triangleBarycenters;
        std::vector<ScalarType> triangleAreas;
        std::vector<CoordType> triangleNormals;
        
        for (auto it = node->oBegin; it != node->oEnd; ++it)
        {
            const FaceType *face = *it;
            if (!face || face->IsD())
                continue;

            // Get triangle vertices
            CoordType v0 = face->cP(0);
            CoordType v1 = face->cP(1);
            CoordType v2 = face->cP(2);

            // Compute triangle area and normal
            CoordType edge01 = v1 - v0;
            CoordType edge02 = v2 - v0;
            CoordType crossProd = edge01 ^ edge02;
            ScalarType area2 = crossProd.Norm();
            ScalarType area = area2 * 0.5;

            if (area < 1e-10)
                continue;

            // Area-weighted normal (unnormalized)
            CoordType N = crossProd * 0.5;
            CoordType barycenter = (v0 + v1 + v2) / 3.0;
            
            // Accumulate
            data.N += N;
            data.area += area;
            data.areaP += barycenter * area;
            
            triangleBarycenters.push_back(barycenter);
            triangleAreas.push_back(area);
            triangleNormals.push_back(N);
        }
        
        // Compute average position
        if (data.area > 1e-10)
        {
            data.averageP = data.areaP / data.area;
        }
        
        // Compute higher-order moments
        // NOTE: For individual triangles at the centroid, Nij = 0
        // The second-order (Nijk) terms are still needed though
#if TAYLOR_SERIES_ORDER >= 1
        data.NijDiag = CoordType(0, 0, 0);
        data.Nxy_Nyx = 0;
        data.Nyz_Nzy = 0;
        data.Nzx_Nxz = 0;
        
#if TAYLOR_SERIES_ORDER >= 2
        // For Nijk, we still need to compute these
        // Using simplified centroid-based approximation
        data.NijkDiag = CoordType(0, 0, 0);
        data.sumPermuteNxyz = 0;
        data.N2xxy_Nyxx = 0;
        data.N2xxz_Nzxx = 0;
        data.N2yyz_Nzyy = 0;
        data.N2yyx_Nxyy = 0;
        data.N2zzx_Nxzz = 0;
        data.N2zzy_Nyzz = 0;
#endif
#endif
        
        // Compute max distance squared from all triangle vertices
        ScalarType maxDist2 = 0.0;
        for (auto it = node->oBegin; it != node->oEnd; ++it)
        {
            const FaceType *face = *it;
            if (!face || face->IsD()) continue;
            
            for (int vi = 0; vi < 3; ++vi)
            {
                ScalarType d = (face->cP(vi) - data.averageP).SquaredNorm();
                maxDist2 = std::max(maxDist2, d);
            }
        }
        data.maxPDist2 = maxDist2;
    }
    else
    {
        // Internal node - recurse to children first
        for (int i = 0; i < 2; ++i)
        {
            if (node->children[i])
                precomputeNode(node->children[i]);
        }
        precomputeInternal(node);
    }
}

// Evaluate winding number using Taylor series approximation
template<typename MeshType>
typename FastWindingNumber<MeshType>::ScalarType 
FastWindingNumber<MeshType>::evaluateApproximation(const CoordType &queryPoint, const NodeType *node) const
{
    if (!node)
        return 0.0;

    const BoxDataType &data = node->auxData;

    // q = averageP - queryPoint (vector from query to center)
    CoordType q = data.averageP - queryPoint;
    ScalarType q2 = q.SquaredNorm();
    
    if (q2 < 1e-10)
        return 0.0;

    ScalarType qNorm = sqrt(q2);
    ScalarType q_inv = 1.0 / qNorm;
    ScalarType q2_inv = 1.0 / q2;
    ScalarType q3_inv = q2_inv * q_inv;
    
    // Normalize q direction
    CoordType qNormalized = q * q_inv;

    // Omega_0 (monopole term): -q·N / |q|^3
    ScalarType omega = -q3_inv * (qNormalized * data.N);

#if TAYLOR_SERIES_ORDER >= 1
    if (myOrder >= 1)
    {
        ScalarType qx = qNormalized[0], qy = qNormalized[1], qz = qNormalized[2];
        ScalarType qx2 = qx * qx, qy2 = qy * qy, qz2 = qz * qz;
        
        ScalarType q5_inv = q3_inv * q2_inv;
        
        // Omega_1 = (1/|q|^5) * [sum(Nii) - 3*(sum(Nii*qi^2) + Nxy*qx*qy + Nyz*qy*qz + Nzx*qz*qx)]
        ScalarType trace = data.NijDiag[0] + data.NijDiag[1] + data.NijDiag[2];
        ScalarType quadratic_term = data.NijDiag[0] * qx2 + data.NijDiag[1] * qy2 + data.NijDiag[2] * qz2 +
                          data.Nxy_Nyx * qx * qy + data.Nyz_Nzy * qy * qz + data.Nzx_Nxz * qz * qx;
        
        ScalarType omega1 = q5_inv * (trace - 3.0 * quadratic_term);
        omega += omega1;

#if TAYLOR_SERIES_ORDER >= 2
        if (myOrder >= 2)
        {
            ScalarType qx3 = qx2 * qx, qy3 = qy2 * qy, qz3 = qz2 * qz;
            ScalarType q7_inv = q5_inv * q2_inv;
            
            // temp0 and temp1 arrays from UT_SolidAngle.cpp
            ScalarType temp0_x = data.N2yyx_Nxyy + data.N2zzx_Nxzz;
            ScalarType temp0_y = data.N2zzy_Nyzz + data.N2xxy_Nyxx;
            ScalarType temp0_z = data.N2xxz_Nzxx + data.N2yyz_Nzyy;
            
            ScalarType temp1_x = qy * data.N2xxy_Nyxx + qz * data.N2xxz_Nzxx;
            ScalarType temp1_y = qz * data.N2yyz_Nzyy + qx * data.N2yyx_Nxyy;
            ScalarType temp1_z = qx * data.N2zzx_Nxzz + qy * data.N2zzy_Nyzz;
            
            ScalarType linear_term = qx * (3.0 * data.NijkDiag[0] + temp0_x) +
                           qy * (3.0 * data.NijkDiag[1] + temp0_y) +
                           qz * (3.0 * data.NijkDiag[2] + temp0_z);
            
            ScalarType cubic_term = data.NijkDiag[0] * qx3 + data.NijkDiag[1] * qy3 + data.NijkDiag[2] * qz3 +
                          data.sumPermuteNxyz * qx * qy * qz +
                          qx2 * temp1_x + qy2 * temp1_y + qz2 * temp1_z;
            
            ScalarType omega2 = q7_inv * (1.5 * linear_term - 7.5 * cubic_term);
            omega += omega2;
        }
#endif
    }
#endif

    return omega;
}

// Traverse the tree and accumulate solid angles
template<typename MeshType>
typename FastWindingNumber<MeshType>::ScalarType 
FastWindingNumber<MeshType>::traverseTree(const CoordType &queryPoint, const NodeType *node, const ScalarType accuracyScale2) const
{
    if (!node)
        return 0.0;

    const BoxDataType &data = node->auxData;

    // Check if we can use approximation
    CoordType r = data.averageP - queryPoint;
    ScalarType r2 = r.SquaredNorm();

    // Accuracy criterion: alpha^2 = (R^2) / (r^2) < threshold^2
    // where R is the bounding box size and r is the distance to query point
    bool useApproximation = (data.maxPDist2 / r2) < accuracyScale2;

    if (useApproximation && !node->IsLeaf())
    {
        // Use Taylor series approximation
        return evaluateApproximation(queryPoint, node);
    }
    else if (node->IsLeaf())
    {
        // Leaf node - compute exact solid angles for all triangles
        ScalarType sum = 0.0;
        for (auto it = node->oBegin; it != node->oEnd; ++it)
        {
            sum += computeTriangleSolidAngle(queryPoint, *it);
        }
        return sum;
    }
    else
    {
        // Recurse to children
        ScalarType sum = 0.0;
        for (int i = 0; i < 2; ++i)
        {
            if (node->children[i])
                sum += traverseTree(queryPoint, node->children[i], accuracyScale2);
        }
        return sum;
    }
}

// Initialize the fast winding number structure
template<typename MeshType>
void FastWindingNumber<MeshType>::init(const MeshType &mesh, const int order)
{
    clear();

    pMesh = &mesh;
    myOrder = order;
    nTriangles = mesh.FN();

    // Build the AABB tree
    std::vector<FacePointer> facePtrs;
    facePtrs.reserve(nTriangles);
    for (size_t i = 0; i < mesh.face.size(); ++i)
    {
        if (!mesh.face[i].IsD())
            facePtrs.push_back(const_cast<FacePointer>(&mesh.face[i]));
    }

    // Define functors required by AABBBinaryTree::Set
    struct GetPointerFunctor {
        FacePointer operator()(FacePointer fp) const { return fp; }
    };
    
    struct GetBoxFunctor {
        void operator()(const FaceType &f, Box3<ScalarType> &box) const {
            box.SetNull();
            box.Add(f.cP(0));
            box.Add(f.cP(1));
            box.Add(f.cP(2));
        }
    };
    
    struct GetBarycenterFunctor {
        void operator()(const FaceType &f, CoordType &barycenter) const {
            barycenter = (f.cP(0) + f.cP(1) + f.cP(2)) / 3.0;
        }
    };

    GetPointerFunctor getPtr;
    GetBoxFunctor getBox;
    GetBarycenterFunctor getBarycenter;

    const unsigned int maxObjectsPerLeaf = 1;
    const ScalarType leafBoxMaxVolume = (ScalarType)0;
    const bool useVariance = true;

    tree.Set(facePtrs.begin(), facePtrs.end(), getPtr, getBox, getBarycenter, 
             maxObjectsPerLeaf, leafBoxMaxVolume, useVariance);

    // Precompute all node data
    precomputeNode(tree.pRoot);
}

// Compute solid angle (winding number contribution) for a query point
template<typename MeshType>
typename FastWindingNumber<MeshType>::ScalarType 
FastWindingNumber<MeshType>::computeSolidAngle(const CoordType &queryPoint, const ScalarType accuracyScale) const
{
    ScalarType accuracyScale2 = accuracyScale * accuracyScale;
    ScalarType windingNumber = traverseTree(queryPoint, tree.pRoot, accuracyScale2);
    return windingNumber;
}

} // namespace tri
} // namespace vcg

#endif // __VCG_TRI_FAST_WINDING_NUMBER
