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
#ifndef __VCGLIB_CURVE_ON_SURF_H
#define __VCGLIB_CURVE_ON_SURF_H

#include<vcg/complex/complex.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/math/histogram.h>
#include<vcg/space/distance3.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <wrap/io_trimesh/export_ply.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh
/// \brief A class for managing curves on a 2-manifold (Curve on Manifold - CoM).
/**
 * This class is used to project, simplify, smooth, and snap polylines (represented as edge meshes) 
 * over a triangulated surface (the "base mesh").
 * 
 * \par Overview
 * The CoM class provides tools to:
 * - Project polylines onto a surface
 * - Snap polyline vertices to mesh vertices or edges
 * - Refine polylines to follow surface features
 * - Simplify polylines while maintaining surface fidelity
 * - Split the base mesh along polylines for mesh cutting operations
 * 
 * \par Terminology
 * - **Base mesh**: The triangulated surface mesh (stored in `base`)
 * - **Polyline/Curve**: An edge mesh passed to the various methods (typically named `poly` in parameters)
 * - **Snapping**: The process of aligning polyline vertices to mesh vertices or edges using barycentric thresholds
 * 
 * \par Requirements
 * The base mesh should be:
 * - 2-manifold
 * - Have reasonable triangle quality
 * - Have up-to-date topology (FaceFace adjacency)
 * - Have bounding box correctly set
 * 
 * \par Usage Pattern
 * 1. Initialize the class with a base mesh: `CoM<MeshType> com(baseMesh);`
 * 2. Call `Init()` to build the spatial acceleration structures
 * 3. Pass polylines (as edge meshes) to various methods for processing
 * 4. Adjust parameters via the `par` member for fine control
 * 
 * \par Implementation Notes
 * - The class uses barycentric coordinates to determine if a point of the polyline should snap to vertices or edges
 * - All spatial queries use a uniform grid for acceleration
 * - Many operations are iterative and may require multiple passes
 * 
 * \note There is some naming inconsistency: methods use both "Curve" and "Polyline" 
 *       interchangeably to refer to the edge mesh being processed.
 * 
 */

template <class MeshType>
class CoM
{
public:
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType      CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef Box3<ScalarType>                  Box3Type;
  typedef Segment3<ScalarType>              Segment3Type;  
  typedef typename vcg::GridStaticPtr<FaceType, ScalarType> MeshGrid;  
  typedef typename vcg::GridStaticPtr<EdgeType, ScalarType> EdgeGrid;
  typedef typename face::Pos<FaceType> PosType;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  
  /**
   * \brief Parameter class controlling the behavior of CoM algorithms
   * 
   * This class contains all the thresholds and tolerances used by the various
   * curve-on-manifold operations. Default values are computed relative to the
   * bounding box diagonal of the base mesh.
   */
  class Param 
  {
  public:
    
    ScalarType surfDistThr;        ///< Max distance between surface and curve; used in simplify and refine
    ScalarType minRefEdgeLen;      ///< Minimal admitted edge length (used in refine: never make edges shorter than this) 
    ScalarType maxSimpEdgeLen;     ///< Maximal admitted edge length (used in simplify: never make edges longer than this) 
    ScalarType maxMoveDelta;       ///< The maximum movement admitted during MoveAndProject (before projection) 
    ScalarType maxSnapThr;         ///< The maximum distance allowed when snapping a polyline vertex onto a mesh vertex (currently unused)
    ScalarType gridBailout;        ///< The maximum distance bailout used in grid-based spatial queries
    ScalarType barycentricSnapThr; ///< Threshold for snapping barycentric coords to 0 or 1 (controls vertex/edge snapping)
    
    /// Constructor with default parameter initialization based on mesh size
    Param(MeshType &m) { SetDefault(m);}
    
    /// Set all parameters to reasonable defaults based on the mesh bounding box
    void SetDefault(MeshType &m)
    {
      surfDistThr        = m.bbox.Diag()/1000.0;
      minRefEdgeLen      = m.bbox.Diag()/16000.0;
      maxSimpEdgeLen     = m.bbox.Diag()/100.0;
      maxMoveDelta       = m.bbox.Diag()/100.0;
      maxSnapThr         = m.bbox.Diag()/1000.0;
      gridBailout        = m.bbox.Diag()/20.0;
      barycentricSnapThr = 0.05;
    }
    
    /// Print current parameter values to stdout
    void Dump() const
    {
      printf("surfDistThr    = %6.4f\n",surfDistThr   );
      printf("minRefEdgeLen  = %6.4f\n",minRefEdgeLen    );
      printf("maxSimpEdgeLen = %6.4f\n",maxSimpEdgeLen    );
      printf("maxMoveDelta   = %6.4f\n",maxMoveDelta);
    }
  };
  
  
  
  // ============================================================================
  // Data Members
  // ============================================================================
  
  MeshType &base;       ///< Reference to the base triangulated surface mesh
  MeshGrid uniformGrid; ///< Spatial acceleration structure for closest point queries
  Param par;            ///< Parameters controlling algorithm behavior
  
  /// Constructor: initializes the CoM with a base mesh
  CoM(MeshType &_m) :base(_m),par(_m){}
 
  // ============================================================================
  // Spatial Query Methods
  // ============================================================================
  
  /**
   * \brief Get the closest face to a query point
   * \param p The query point
   * \return Pointer to the closest face, or nullptr if none found within gridBailout distance
   */
  FaceType *GetClosestFace(const CoordType &p)
  {
    ScalarType closestDist;
    CoordType closestP;
    return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP);
  }
  
  /**
   * \brief Get the closest face and its barycentric coordinates
   * \param p The query point
   * \param ip Output: barycentric coordinates of the closest point on the returned face
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip)
    {
      ScalarType closestDist;
      CoordType closestP,closestN;
      return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP,closestN,ip);
    }

  /**
   * \brief Get the closest face, barycentric coordinates, and normal
   * \param p The query point
   * \param ip Output: barycentric coordinates of the closest point on the returned face
   * \param in Output: normal at the closest point
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip, CoordType &in)
    {
      ScalarType closestDist;
      CoordType closestP;
      return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP,in,ip);
    }

  /**
   * \brief Get the closest face and the closest point on it
   * \param p The query point
   * \param closestP Output: the 3D coordinates of the closest point on the surface
   * \return Pointer to the closest face
   */
  FaceType *GetClosestFacePoint(const CoordType &p, CoordType &closestP)
  {
    ScalarType closestDist;
    return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP);
  }
  
  // ============================================================================
  // Barycentric Coordinate and Snapping Methods
  // ============================================================================
  
  /** 
   * \brief Test if a barycentric coordinate is well snapped 
   * \param ip The barycentric coordinate to test (must sum to 1.0)
   * \return true if no snapping is needed (all coords are either 0, 1, or far from boundaries)
   * 
   * A barycentric coordinate is "well snapped" if each component is either:
   * - Exactly 0.0 or 1.0, OR
   * - Far enough from 0 and 1 (outside the barycentricSnapThr threshold)
   * 
   * This indicates the point doesn't need further snapping adjustment.
   */
  bool IsWellSnapped(const CoordType &ip)
  {
      for(int i=0;i<3;++i)
          if( (ip[i]< par.barycentricSnapThr         && ip[i]!= 0.0) ||
              (ip[i]> (1.0 - par.barycentricSnapThr) && ip[i]!= 1.0))
              return false;
      assert(ip[0]+ip[1]+ip[2] == 1.0);
      return true;
  }
  
  /**
   * \brief Check if a barycentric coordinate is snapped to an edge
   * \param ip The barycentric coordinate (must be snapped)
   * \param ei Output: index (0,1,2) of the edge opposite to the zero coordinate
   * \return true if snapped to an edge (exactly one coordinate is 0, other two are positive)
   * 
   * Edge snapping means the point lies on one of the three edges of the triangle.
   * Edge i is the edge from vertex i to vertex (i+1)%3, opposite to vertex (i+2)%3.
   */
  bool IsSnappedEdge(CoordType &ip, int &ei)
  {
    for(int i=0;i<3;++i)
      if(ip[i]>0.0 && ip[(i+1)%3]>0.0 && ip[(i+2)%3]==0.0 ) {
        ei=i;
        return true; 
      }
    ei=-1;
    return false;
  }

  /**
   * \brief Check if a barycentric coordinate is snapped to a vertex
   * \param ip The barycentric coordinate (must be snapped)
   * \param vi Output: index (0,1,2) of the vertex with coordinate == 1.0
   * \return true if snapped to a vertex (one coordinate is 1.0, others are 0.0)
   */
  bool IsSnappedVertex(CoordType &ip, int &vi)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) {
        vi=i;
        return true; 
      }
    vi=-1;
    return false;
  }

  /**
   * \brief Find the vertex pointer for a vertex-snapped barycentric coordinate
   * \param fp The face containing the point
   * \param ip The barycentric coordinate (should be snapped to a vertex)
   * \return Pointer to the snapped vertex, or nullptr if not vertex-snapped
   */
  VertexPointer FindVertexSnap(FacePointer fp, CoordType &ip)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) return fp->V(i);
    return 0;
  }
  
  // ============================================================================
  // Polyline Tagging and Mesh Cutting Methods
  // ============================================================================
  
  /**
   * \brief Tag face edges of the base mesh that coincide with polyline edges
   * \param poly The polyline (as edge mesh) to use for tagging
   * \param markFlag If true, clears all FaceEdgeS flags before tagging (default: true)
   * \return true if ALL edges of the polyline are fully snapped onto mesh edges
   * 
   * This function marks edges in the base mesh (using FaceEdgeS flag) where they 
   * coincide with edges from the polyline. The polyline edges must be snapped to 
   * mesh vertices at both endpoints for this to work.
   * This function requires VertexFace and FaceFace adjacency.
   * 
   * \note This is typically used as a preparation step before cutting the mesh
   *       along the polyline using CutMeshAlongCrease or similar functions.
   *      
   * 
   * \warning Returns false if any polyline edge is not properly snapped, or if
   *          the snapped vertices don't form an edge in the base mesh.
   * 
   * \sa SplitMeshWithPolyline, CutMeshAlongCrease
   */
    
bool TagFaceEdgeSelWithPolyLine(MeshType &poly,bool markFlag=true)
{
	if (markFlag)
		tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(base);

	tri::UpdateTopology<MeshType>::VertexFace(base);
	tri::UpdateTopology<MeshType>::FaceFace(base);

	for(EdgeIterator ei=poly.edge.begin(); ei!=poly.edge.end();++ei)
	{
		CoordType ip0,ip1;
		FaceType *f0 = GetClosestFaceIP(ei->cP(0),ip0);
		FaceType *f1 = GetClosestFaceIP(ei->cP(1),ip1);

		if(BarycentricSnap(ip0) && BarycentricSnap(ip1))
		{
			VertexPointer v0 = FindVertexSnap(f0,ip0);
			VertexPointer v1 = FindVertexSnap(f1,ip1);

			if(v0==0 || v1==0)
				return false;
			if(v0==v1)
				return false;

			FacePointer ff0,ff1;
			int e0,e1;
			bool ret=face::FindSharedFaces<FaceType>(v0,v1,ff0,ff1,e0,e1);
			if(ret)
			{
				assert(ret);
				assert(ff0->V(e0)==v0 || ff0->V(e0)==v1);
				ff0->SetFaceEdgeS(e0);
				ff1->SetFaceEdgeS(e1);
			} else {
				return false;
			}
		}
		else {
			return false;
		}
	}
	return true;
}

  /**
   * \brief Find the minimum distance from a sample point to the polyline
   * \param samplePnt The point to measure distance from
   * \param edgeGrid Spatial acceleration grid for the polyline edges
   * \param poly The polyline (as edge mesh)
   * \param closestPoint Output: the closest point on the polyline
   * \return The minimum distance from samplePnt to the polyline
   */
  ScalarType MinDistOnEdge(CoordType samplePnt, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
      ScalarType polyDist;
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,par.gridBailout,polyDist,closestPoint);        
      return polyDist;    
  }
  
  /**
   * \brief Find the closest point on a mesh edge to the polyline (static version)
   * \param v0 First vertex of the mesh edge
   * \param v1 Second vertex of the mesh edge
   * \param edgeGrid Spatial acceleration grid for the polyline edges
   * \param poly The polyline (as edge mesh)
   * \param closestPoint Output: the point on the edge [v0,v1] closest to the polyline
   * \return The minimum distance from the edge to the polyline
   * 
   * This samples the edge [v0,v1] uniformly and finds which sample is closest to the polyline.
   */
  static ScalarType MinDistOnEdge(VertexType *v0,VertexType *v1, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
    ScalarType minPolyDist = std::numeric_limits<ScalarType>::max();
    const ScalarType sampleNum = 50;
    const ScalarType maxDist = poly.bbox.Diag()/10.0;
    for(ScalarType k = 0;k<sampleNum+1;++k)
    {
      ScalarType polyDist;
      CoordType closestPPoly;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,maxDist,polyDist,closestPPoly);        
      
      if(polyDist < minPolyDist)
      {
        minPolyDist = polyDist;
        closestPoint = samplePnt;
//        closestPoint = closestPPoly;
      }
    }
    return minPolyDist;    
  }
  
  // ============================================================================
  // Attribute Extraction and Comparison (for Seam Processing)
  // ============================================================================
  
  /**
   * \brief Extract vertex attributes for seam processing
   * \param srcMesh Source mesh (unused but required by interface)
   * \param f The face containing the vertex
   * \param whichWedge Which vertex (0,1,2) of the face to extract
   * \param dstMesh Destination mesh (unused but required by interface)
   * \param v Output: vertex with copied attributes
   * 
   * This is a callback function used by the attribute_seam system.
   * It copies all per-vertex properties and uses the face color.
   * 
   * \note This is used when splitting the mesh along seams/polylines.
   */
  static inline void ExtractVertex(const MeshType & srcMesh, const FaceType & f, int whichWedge, const MeshType & dstMesh, VertexType & v)
  {
      (void)srcMesh;
      (void)dstMesh;
      // This is done to preserve every single perVertex property
      // perVextex Texture Coordinate is instead obtained from perWedge one.
      v.ImportData(*f.cV(whichWedge));
      v.C() = f.cC();
  }
  
  /**
   * \brief Compare two vertices for seam compatibility
   * \param m The mesh (unused but required by interface)
   * \param vA First vertex
   * \param vB Second vertex
   * \return true if vertices are compatible across a seam
   * 
   * This callback is used by the attribute_seam system to determine if two
   * vertices can be considered the same across a seam boundary.
   * Current implementation: Red and Blue colored vertices are considered incompatible.
   * 
   * \note This is part of the mesh cutting/seam processing infrastructure.
   */
  static inline bool CompareVertex(const MeshType & m, const VertexType & vA, const VertexType & vB)
  {
      (void)m;
      
      if(vA.C() == Color4b(Color4b::Red) && vB.C() == Color4b(Color4b::Blue) ) return false;
      if(vA.C() == Color4b(Color4b::Blue) && vB.C() == Color4b(Color4b::Red) ) return false;
      return true;      
  }
  
  // ============================================================================
  // Utility Functions
  // ============================================================================
  
  /**
   * \brief Compute quality-weighted linear interpolation between two vertices
   * \param v0 First vertex
   * \param v1 Second vertex
   * \return Interpolated position weighted by inverse quality values
   * 
   * Points with higher quality (larger absolute value) contribute less to the result.
   * This is useful for adaptive refinement based on error metrics stored in quality.
   */
  static CoordType QLerp(VertexType *v0, VertexType *v1)
  {
    
    ScalarType qSum = fabs(v0->Q())+fabs(v1->Q());      
    ScalarType w0 = (qSum - fabs(v0->Q()))/qSum;
    ScalarType w1 = (qSum - fabs(v1->Q()))/qSum;      
    return v0->P()*w0 + v1->P()*w1;      
  }
  
  
  /**
   * @brief SnapPolyline snaps the vertexes of a polyline onto the base mesh
   * @param poly
   * @param newVertVec the vector of the indexes of the snapped vertices
   * @return true if it has modified the polyline
   * 
   * Polyline vertices can be snapped either on vertexes or on edges. 
   * Usually the only points that we should allow to not be snapped are the endpoints and non manifold points.
   * Vertexes are colored according to their snapping state 
   * 
   */  
    
  bool SnapPolyline(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int vertSnapCnt=0;
    int edgeSnapCnt=0;
    int borderCnt=0,midCnt=0,nonmanifCnt=0;
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(BarycentricSnap(ip))
      {
        if(ip[0]>0 && ip[1]>0) { vi->P() = f->P(0)*ip[0]+f->P(1)*ip[1]; edgeSnapCnt++; assert(ip[2]==0); vi->C()=Color4b::White;}
        if(ip[0]>0 && ip[2]>0) { vi->P() = f->P(0)*ip[0]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[1]==0); vi->C()=Color4b::White;}
        if(ip[1]>0 && ip[2]>0) { vi->P() = f->P(1)*ip[1]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[0]==0); vi->C()=Color4b::White;}
        
        if(ip[0]==1.0) { vi->P() = f->P(0); vertSnapCnt++; assert(ip[1]==0 && ip[2]==0); vi->C()=Color4b::Black;  }
        if(ip[1]==1.0) { vi->P() = f->P(1); vertSnapCnt++; assert(ip[0]==0 && ip[2]==0); vi->C()=Color4b::Black;}
        if(ip[2]==1.0) { vi->P() = f->P(2); vertSnapCnt++; assert(ip[0]==0 && ip[1]==0); vi->C()=Color4b::Black;}
      }
      else
      {
        int deg = edge::VEDegree<EdgeType>(&*vi);
        if (deg > 2) { nonmanifCnt++; vi->C()=Color4b::Magenta; }
        if (deg < 2) { borderCnt++;   vi->C()=Color4b::Green;} 
        if (deg== 2) { midCnt++;      vi->C()=Color4b::Blue;} 
      }
    }
    printf("SnapPolyline %i vertices:  snapped %i onto vert and %i onto edges %i nonmanif, %i border, %i mid\n",
           poly.vn, vertSnapCnt, edgeSnapCnt, nonmanifCnt,borderCnt,midCnt); fflush(stdout);
    int dupCnt=tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    if(dupCnt) printf("SnapPolyline: Removed %i Duplicated vertices\n",dupCnt);
    
    return vertSnapCnt==0 && edgeSnapCnt==0 && dupCnt==0;
  }
  
   void SelectBoundaryVertex(MeshType &poly)
   {
     tri::UpdateSelection<MeshType>::VertexClear(poly);
     tri::UpdateTopology<MeshType>::VertexEdge(poly);
     ForEachVertex(poly, [&](VertexType &v){
       if(edge::VEDegree<EdgeType>(&v)==1) v.SetS();
     });
   }
  
   void SelectUniformlyDistributed(MeshType &poly, int k)
   {
     tri::TrivialPointerSampler<MeshType> tps;
     ScalarType samplingRadius = tri::Stat<MeshType>::ComputeEdgeLengthSum(poly)/ScalarType(k);
     tri::SurfaceSampling<MeshType, typename tri::TrivialPointerSampler<MeshType> >::EdgeMeshUniform(poly,tps,samplingRadius);     
     for(int i=0;i<tps.sampleVec.size();++i)
       tps.sampleVec[i]->SetS();
   }
   
   
    
  /*
   * Make an edge mesh 1-manifold by splitting all the
   * vertexes that have more than two incident edges
   * 
   * It performs the split in three steps. 
   * - First it collects and counts the vertices to be splitten. 
   * - Then it adds the vertices to the mesh and 
   * - lastly it updates the poly with the newly added vertices. 
   * 
   * singSplitFlag allows to ubersplit each singularity in a number of vertex of the same order of its degree. 
   * This is not really necessary but helps the management of sharp turns in the poly mesh.
   * \todo add corner detection and split.
   */
  
  void DecomposeNonManifoldPolyline(MeshType &poly, bool singSplitFlag = true)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    std::vector<int> degreeVec(poly.vn, 0);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int neededVert=0;
    int delta;
    if(singSplitFlag) delta = 1;
                 else delta = 2;
      
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      std::vector<EdgeType *> starVec;
      edge::VEStarVE(&*vi,starVec);
      degreeVec[tri::Index(poly, *vi)] = starVec.size();
      if(starVec.size()>2)
        neededVert += starVec.size()-delta;
    }
    printf("DecomposeNonManifold Adding %i vert to a polyline of %i vert\n",neededVert,poly.vn);
    VertexIterator firstVi = tri::Allocator<MeshType>::AddVertices(poly,neededVert);
    
    for(size_t i=0;i<degreeVec.size();++i)
    {
      if(degreeVec[i]>2)
      {
        std::vector<EdgeType *> edgeStarVec;
        edge::VEStarVE(&(poly.vert[i]),edgeStarVec);
        assert(edgeStarVec.size() == degreeVec[i]);
        for(size_t j=delta;j<edgeStarVec.size();++j)
        {
          EdgeType *ep = edgeStarVec[j];
          int ind; // index of the vertex to be changed
          if(tri::Index(poly,ep->V(0)) == i) ind = 0;
              else ind = 1;
  
          ep->V(ind) = &*firstVi;
          ep->V(ind)->P() = poly.vert[i].P();
          ep->V(ind)->N() = poly.vert[i].N();
          ++firstVi;
        }
      }
    }
    assert(firstVi == poly.vert.end());
  }
  
  // ============================================================================
  // Mesh Splitting and Refinement for Polyline Integration
  // ============================================================================
  
  /**
   * \brief Split the base mesh to make it conforming with the polyline
   * \param poly The polyline to integrate into the base mesh
   * 
   * Note that you have to call 
   * RefineCurveByBaseMesh before this function to ensure that the polyline 
   * is sufficiently sampled and follows the surface features.
   * This is a complex two-phase algorithm:
   * 
   * **Phase 1: Vertex Insertion**
   * - Finds all polyline vertices that are NOT snapped to mesh vertices/edges
   * - Inserts them into the base mesh using 1-to-3 face splits
   * - This creates new vertices in the mesh at polyline vertex positions
   * 
   * **Phase 2: Edge Splitting** (iterative)
   * - Identifies mesh edges that are crossed by polyline edges
   * - Splits those edges at the intersection points
   * - Repeats until no more edges need splitting
   * 
   * After completion, the polyline edges should coincide with base mesh edges,
   * allowing subsequent operations like mesh cutting or constrained triangulation.
   * 
   * \note This modifies the base mesh topology! 
   * \note Calls Init() internally to rebuild spatial structures
   * \note Calls SnapPolyline() to update polyline after mesh modifications
   * 
   * \warning This can significantly increase mesh complexity
   * 
   * \sa TagFaceEdgeSelWithPolyLine, SnapPolyline, RefineCurveByBaseMesh
   */
  
  void SplitMeshWithPolyline(MeshType &poly)
  {        
    std::vector< std::pair<int,VertexPointer> > toSplitVec;  // the index of the face to be split and the poly vertex to be used
    
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(!BarycentricSnap(ip))
        toSplitVec.push_back(std::make_pair(tri::Index(base,f),&*vi));
    }
    SimplifyNullEdges(poly);
    printf("SplitMeshWithPolyline found %lu non snapped points\n",toSplitVec.size()); fflush(stdout);

    FaceIterator newFi = tri::Allocator<MeshType>::AddFaces(base,toSplitVec.size()*2);
    VertexIterator newVi = tri::Allocator<MeshType>::AddVertices(base,toSplitVec.size());    
    tri::UpdateColor<MeshType>::PerVertexConstant(base,Color4b::White);
    
    for(size_t i =0; i<toSplitVec.size();++i)
    {
        newVi->P() = toSplitVec[i].second->P(); 
        newVi->C()=Color4b::Green;      
        face::TriSplit(&base.face[toSplitVec[i].first],&*(newFi++),&*(newFi++),&*(newVi++));
    }
    Init(); //  need to reset everthing
    SnapPolyline(poly);
    
    // Second loop to perform the face-face Edge split **********************
    // This loop must be iterated multiple times becouse it can happen that more than one polyline vertices falls on the same edge.
    // So multiple splits must be done.
    std::map<std::pair<CoordType,CoordType>, VertexPointer> edgeToSplitMap;
    do
    {
        edgeToSplitMap.clear();
        for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
        {
            CoordType ip;
            FaceType *f = GetClosestFaceIP(vi->cP(),ip);
            if(!BarycentricSnap(ip)) { assert(0); }            
            for(int i=0;i<3;++i)
            {
                if((ip[i      ]>0 && ip[i      ]<1.0) &&
                    (ip[(i+1)%3]>0 && ip[(i+1)%3]<1.0) &&
                    ip[(i+2)%3]==0 )
                {
                    CoordType p0=f->P0(i);
                    CoordType p1=f->P1(i);
                    if (p0>p1) std::swap(p0,p1);
                    if(edgeToSplitMap[std::make_pair(p0,p1)])
                        printf("Found an already used Edge %lu - %lu vert %lu!!!\n", tri::Index(base,f->V0(i)),tri::Index(base,f->V1(i)),tri::Index(poly,&*vi));
                    edgeToSplitMap[std::make_pair(p0,p1)]=&*vi;
                }         
            }
        }
        printf("SplitMeshWithPolyline: Created a map of %lu edges to be split\n",edgeToSplitMap.size());
        EdgePointPred ePred(edgeToSplitMap);
        EdgePointSplit eSplit(edgeToSplitMap);
        tri::UpdateTopology<MeshType>::FaceFace(base);
        tri::RefineE(base,eSplit,ePred);     
        Init(); //  need to reset everthing
    } while(edgeToSplitMap.size()>0); // while there are edges to be split
 }
    
  // ============================================================================
  // Initialization
  // ============================================================================
  
  /**
   * \brief Initialize the CoM data structures for processing
   * 
   * This must be called after construction and whenever the base mesh is modified.
   * It performs:
   * - Face normal computation
   * - Face-Face topology update
   * - Spatial acceleration structure (uniform grid) construction
   * 
   * \note Call this before using any polyline processing methods.
   * \warning If the base mesh topology changes, call Init() again.
   */
  void Init()
  {
    // Construction of the uniform grid
    UpdateNormal<MeshType>::PerFaceNormalized(base);
    UpdateTopology<MeshType>::FaceFace(base);    
    uniformGrid.Set(base.face.begin(), base.face.end());    
  }
  
  // ============================================================================
  // Simplification Methods
  // ============================================================================
  
  /**
   * \brief Remove duplicate/zero-length edges from a polyline
   * \param poly The polyline to simplify
   * 
   * Removes vertices that have collapsed to the same position.
   */
  void SimplifyNullEdges(MeshType &poly)
  {
      int cnt=tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(cnt)
          printf("SimplifyNullEdges: Removed %i Duplicated vertices\n",cnt);
  }
  
  void SimplifyMidEdge(MeshType &poly)
  {
   int startVn;
   int midEdgeCollapseCnt=0;
   tri::Allocator<MeshType>::CompactEveryVector(poly); 
   do
   {
    startVn = poly.vn;
    for(int ei =0; ei<poly.en; ++ei)
    {
      VertexType *v0=poly.edge[ei].V(0);
      VertexType *v1=poly.edge[ei].V(1);
      CoordType ip0,ip1;    
      FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
      FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
      
      bool snap0=BarycentricSnap(ip0);
      bool snap1=BarycentricSnap(ip1);
      int e0i,e1i;
      bool e0 = IsSnappedEdge(ip0,e0i);
      bool e1 = IsSnappedEdge(ip1,e1i);
      if(e0 && e1)
        if( (          f0 == f1           &&          e0i == e1i) || 
            (          f0 == f1->FFp(e1i) &&          e0i == f1->FFi(e1i)) || 
            (f0->FFp(e0i) == f1           && f0->FFi(e0i) == e1i) || 
            (f0->FFp(e0i) == f1->FFp(e1i) && f0->FFi(e0i) == f1->FFi(e1i)) ) 
        {
          CoordType newp = (v0->P()+v1->P())/2.0;
          v0->P()=newp;
          v1->P()=newp;
          midEdgeCollapseCnt++;
        }
    }
    tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
//    printf("SimplifyMidEdge %5i -> %5i %i mid %i ve \n",startVn,poly.vn,midEdgeCollapseCnt);
   } while(startVn>poly.vn);
  } 
  
  /**
   * @brief SimplifyMidFace remove all the vertices that in the mid of a face 
   * and between two of the points snapped onto the edges of the same face
   * @param poly
   * 
   * It assumes that the mesh has been snapped and refined by the BaseMesh
   * 
   */
  void SimplifyMidFace(MeshType &poly)
  {
   int startVn= poly.vn;;
   int midFaceCollapseCnt=0;
   int vertexEdgeCollapseCnt=0;
   int curVn;
   do
   {
    tri::Allocator<MeshType>::CompactEveryVector(poly); 
    curVn = poly.vn;
    UpdateTopology<MeshType>::VertexEdge(poly);
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if( (starVecVp.size()==2) )
      {
        CoordType ipP, ipN, ipI; 
        FacePointer fpP = GetClosestFaceIP(starVecVp[0]->P(),ipP);
        FacePointer fpN = GetClosestFaceIP(starVecVp[1]->P(),ipN);
        FacePointer fpI = GetClosestFaceIP(poly.vert[i].P(), ipI);
        
        bool snapP = (BarycentricSnap(ipP));
        bool snapN = (BarycentricSnap(ipN));
        bool snapI = (BarycentricSnap(ipI));
        VertexPointer vertexSnapP = 0;
        VertexPointer vertexSnapN = 0;
        VertexPointer vertexSnapI = 0;
        for(int j=0;j<3;++j)
        {
          if(ipP[j]==1.0) vertexSnapP=fpP->V(j);
          if(ipN[j]==1.0) vertexSnapN=fpN->V(j);
          if(ipI[j]==1.0) vertexSnapI=fpI->V(j);
        }
        
        bool collapseFlag=false;
        
        if((!snapI && snapP && snapN) ||              // First case a vertex that is not snapped between two snapped vertexes 
           (!snapI && !snapP && fpI==fpP) || // Or a two vertex not snapped but on the same face
           (!snapI && !snapN && fpI==fpN) )
        {
          collapseFlag=true;
          midFaceCollapseCnt++;
        } 
        
        else  // case 2) a vertex snap and edge snap we have to check that the edge do not share the same vertex of the vertex snap
          if(snapI && snapP && snapN && vertexSnapI==0 && (vertexSnapP!=0 || vertexSnapN!=0) )
          {
            for(int j=0;j<3;++j) {
              if(ipI[j]!=0 && (fpI->V(j)==vertexSnapP || fpI->V(j)==vertexSnapN)) {
                collapseFlag=true;                                          
                vertexEdgeCollapseCnt++;
              }
            }
          }            
        
        if(collapseFlag)  
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));
      }
    }  
   } while(curVn>poly.vn);
   printf("SimplifyMidFace %5i -> %5i %i mid %i ve \n",startVn,poly.vn,midFaceCollapseCnt,vertexEdgeCollapseCnt);
  } 
  
  void Simplify(MeshType &poly)
  {
    int startEn = poly.en;
    Distribution<ScalarType> hist;
    for(int i =0; i<poly.en;++i) 
      hist.Add(edge::Length(poly.edge[i]));
        
    UpdateTopology<MeshType>::VertexEdge(poly);
    
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if ((starVecVp.size()==2) && (!poly.vert[i].IsS()))
      {
        ScalarType newSegLen = Distance(starVecVp[0]->P(), starVecVp[1]->P());
        Segment3Type seg(starVecVp[0]->P(),starVecVp[1]->P());
        ScalarType segDist;
        CoordType closestPSeg;
        SegmentPointDistance(seg,poly.vert[i].cP(),closestPSeg,segDist);
        CoordType fp,fn;
        ScalarType maxSurfDist = MaxSegDist(starVecVp[0], starVecVp[1],fp,fn);
        
        if((maxSurfDist < par.surfDistThr) && (newSegLen < par.maxSimpEdgeLen) )
        {
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));          
        }
      }
    }
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
//    printf("Simplify %5i -> %5i (total len %5.2f)\n",startEn,poly.en,hist.Sum());
  }
  
  void EvaluateHausdorffDistance(MeshType &poly, Distribution<ScalarType> &dist)
  {
    dist.Clear();
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    tri::UpdateQuality<MeshType>::VertexConstant(poly,0);
    for(int i =0; i<poly.edge.size();++i)
    {      
      CoordType farthestP, farthestN;      
      ScalarType maxDist = MaxSegDist(poly.edge[i].V(0),poly.edge[i].V(1), farthestP, farthestN, &dist);      
      poly.edge[i].V(0)->Q()+= maxDist;
      poly.edge[i].V(1)->Q()+= maxDist;
    }
    for(int i=0;i<poly.vn;++i)
    {
      ScalarType deg = edge::VEDegree<EdgeType>(&poly.vert[i]);
      poly.vert[i].Q()/=deg;
    }
    tri::UpdateColor<MeshType>::PerVertexQualityRamp(poly,0,dist.Max());    
  }
  

  /**
   * \brief Snap barycentric coordinates to 0 or 1 if within threshold
   * \param ip Input/Output: barycentric coordinates (must sum to 1.0)
   * \return true if the point was snapped to a vertex or edge (at least one coord became 0)
   * 
   * **This is one of the MOST IMPORTANT functions in the class** - it's used throughout!
   * 
   * Given barycentric coordinates of a point in a triangle, this function decides 
   * whether it should be "snapped" to a vertex or edge based on the 
   * `par.barycentricSnapThr` threshold.
   * 
   * **Algorithm:**
   * 1. If any coordinate is within `barycentricSnapThr` of 0, snap it to 0
   * 2. If any coordinate is within `barycentricSnapThr` of 1, snap it to 1
   * 3. Renormalize to ensure sum = 1.0
   * 4. If sum is still not exactly 1.0 (due to floating point), adjust the non-snapped coordinate
   * 
   * **Snapping Cases:**
   * - One coord = 1.0, others = 0 → Snapped to a vertex
   * - One coord = 0, others > 0 → Snapped to an edge
   * - All coords > 0 and < 1 → Interior point, NOT snapped
   * 
   * **Return Value:**
   * - `true`: Point is on a vertex or edge (at least one coordinate is 0)
   * - `false`: Point is in the interior of the triangle
   * 
   * \note This function MODIFIES the input coordinates in-place!
   * \note The threshold `par.barycentricSnapThr` (default 0.05) controls snapping sensitivity
   * 
   * \warning Side effect: modifies ip parameter! Consider renaming to BarycentricSnapInPlace()
   * 
   * \sa IsWellSnapped, IsSnappedVertex, IsSnappedEdge
   */
  bool BarycentricSnap(CoordType &ip)
  {
    for(int i=0;i<3;++i)
    {
      if(ip[i] <= par.barycentricSnapThr) ip[i]=0;
      if(ip[i] >= 1.0-par.barycentricSnapThr) ip[i]=1;
    }
    ScalarType sum = ip[0]+ip[1]+ip[2];
    
    for(int i=0;i<3;++i) 
      if(ip[i]!=1.0) ip[i]/=sum;
    
    sum = ip[0]+ip[1]+ip[2];
    
    if(sum!=1.0){
        for(int i=0;i<3;++i)
            if(ip[i]>0.0 && ip[i]<1.0) // if it is non snapped
                ip[i]=1.0-(ip[(i+1)%3]+ip[(i+2)%3]);
    }
    
    sum = ip[0]+ip[1]+ip[2];     
    assert(sum ==1.0);
    assert(IsWellSnapped(ip));
    if(ip[0]==0 || ip[1]==0 || ip[2]==0) return true;
    return false;
  }
  
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMesh(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;
    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=false,snap1=false; // true if the segment start/end on a edge/vert
    
    Segment3Type seg0; // The two segments to be avoided 
    Segment3Type seg1; // from which the current poly segment can start
    VertexPointer vertexSnap0 = 0;
    VertexPointer vertexSnap1 = 0;
    if(BarycentricSnap(ip0)) { 
      snap0=true; 
      for(int i=0;i<3;++i) {
        if(ip0[i]==1.0) vertexSnap0=f0->V(i);
        if(ip0[i]==0.0) seg0=Segment3Type(f0->P1(i),f0->P2(i)); 
      }        
    } 
    if(BarycentricSnap(ip1)) { 
      snap1=true; 
      for(int i=0;i<3;++i){
        if(ip1[i]==1.0) vertexSnap1=f1->V(i);
        if(ip1[i]==0.0) seg1=Segment3Type(f1->P1(i),f1->P2(i)); 
      }        
    } 
    
    CoordType bestSplitPt(0,0,0);
    ScalarType bestDist = std::numeric_limits<ScalarType>::max();
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      CoordType samplePnt = segPoly.Lerp(k/sampleNum);    
      CoordType ip;
      FaceType *f=GetClosestFaceIP(samplePnt,ip);
//      BarycentricEdgeSnap(ip);
      if(BarycentricSnap(ip))
      {
        VertexPointer vertexSnapI = 0;        
        for(int i=0;i<3;++i)
          if(ip[i]==1.0) vertexSnapI=f->V(i);
        CoordType closestPt = f->P(0)*ip[0]+f->P(1)*ip[1]+f->P(2)*ip[2];
        if(Distance(samplePnt,closestPt) < bestDist )  
        {
          ScalarType dist0=std::numeric_limits<ScalarType>::max();
          ScalarType dist1=std::numeric_limits<ScalarType>::max();
          CoordType closestSegPt;
          if(snap0) SegmentPointDistance(seg0,closestPt,closestSegPt,dist0);
          if(snap1) SegmentPointDistance(seg1,closestPt,closestSegPt,dist1);
          if( (!vertexSnapI && (dist0 > par.surfDistThr/1000 && dist1>par.surfDistThr/1000) ) ||
              ( vertexSnapI!=vertexSnap0 && vertexSnapI!=vertexSnap1)  )
          {
            bestDist = Distance(samplePnt,closestPt);
            bestSplitPt = closestPt;            
          }
        }      
      }
    }
    if(bestDist < par.surfDistThr*100)
    {
      splitPt = bestSplitPt;
      return true;
    }
    
    return false;
  }
  /**
   * @brief SnappedOnSameFace Return true if the two points are snapped to a common face;
   * @param f0
   * @param i0
   * @param f1
   * @param i0
   * @return 
   * 
   * Require FFAdj. se assume that both SNAPPED. Three cases:
   * - Edge Edge - true iff the two edges belongs to a common face. 
   * - Vert Edge - true iff there is one of the two snapped edge faces has the vert as non-edge face;  
   * - Vert Vert 
   * 
   */
  bool SnappedOnSameFace(FacePointer f0, CoordType i0, FacePointer f1, CoordType i1)
  {
   if(f0==f1) return true;
   int e0,e1;
   int v0,v1;
   bool e0Snap = IsSnappedEdge(i0,e0);
   bool e1Snap = IsSnappedEdge(i1,e1);
   bool v0Snap = IsSnappedVertex(i0,v0);
   bool v1Snap = IsSnappedVertex(i1,v1);
   FacePointer f0p=0; int e0p=-1;  // When Edge snap the other face and the index of the snapped edge on the other face
   FacePointer f1p=0; int e1p=-1;
   assert((e0Snap != v0Snap) && (e1Snap != v1Snap));
   // For EdgeSnap compute the 'other' face stuff 
   if(e0Snap){
     f0p = f0->FFp(e0); e0p=f0->FFi(e0); assert(f0p->FFp(e0p)==f0);
   }
   if(e1Snap){
     f1p = f1->FFp(e1); e1p=f1->FFi(e1); assert(f1p->FFp(e1p)==f1);
   }
   
   if(e0Snap && e1Snap) {
    if(f0==f1p || f0p==f1p || f0p==f1 || f0==f1) return true;
   }
   
   if(e0Snap && v1Snap)  {
     assert(v1>=0 && v1<3 && v0==-1 && e1==-1);
     if(f0->V2(e0)  ==f1->V(v1)) return true;
     if(f0p->V2(e0p)==f1->V(v1)) return true;
   }
     
   if(e1Snap && v0Snap)  {
     assert(v0>=0 && v0<3 && v1==-1 && e0==-1);
     if(f1->V2(e1)  ==f0->V(v0)) return true;
     if(f1p->V2(e1p)==f0->V(v0)) return true;
   }
     
   if(v1Snap && v0Snap)  {
     PosType startPos(f0,f0->V(v0));
     PosType curPos=startPos;
     do
     {
       assert(curPos.V()==f0->V(v0));
       if(curPos.VFlip()==f1->V(v1)) return true;
       curPos.FlipE();
       curPos.FlipF();       
     }
     while(curPos!=startPos);   
   }   
   return false;    
  }
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMeshAdapt(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    splitPt=(v0->P()+v1->P())/2.0;
      
    CoordType ip0,ip1,ipm;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    FaceType *fm=GetClosestFaceIP(splitPt,ipm);
    
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    bool snapm=BarycentricSnap(ipm);
    
    splitPt = fm->P(0)*ipm[0]+fm->P(1)*ipm[1]+fm->P(2)*ipm[2];
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;            
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    if(snap1) {
      int e1,v1;
      if (IsSnappedEdge(ip1,e1)) {
        if(f1->FFp(e1) == f0) return false;
      }
      if(IsSnappedVertex(ip1,v1)) {
        for(int i=0;i<3;++i) 
          if(f0->V(i)==f1->V(v1)) return false;
      }
    }
    
    return true;
  }
  
  
  bool TestSplitSegWithMeshAdaptOld(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      splitPt=(v0->P()+v1->P())/2.0;
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;      
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    splitPt=(v0->P()+v1->P())/2.0;
    return true;
  }
  
  // Given a segment find the maximum distance from it to the original surface. 
  // It is used to evaluate the Haustdorff distance of a Segment from the mesh.
  ScalarType MaxSegDist(VertexType *v0, VertexType *v1, CoordType &farthestPointOnSurf, CoordType &farthestN, Distribution<ScalarType> *distanceDistribution=0)
  {
    ScalarType maxSurfDist = 0;
    const ScalarType sampleNum = 10;
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      ScalarType surfDist;
      CoordType closestPSurf;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,samplePnt,par.gridBailout, surfDist, closestPSurf);        
      if(distanceDistribution)
        distanceDistribution->Add(surfDist);
      assert(f);
      if(surfDist > maxSurfDist)
      {
        maxSurfDist = surfDist;
        farthestPointOnSurf = closestPSurf;
        farthestN = f->N();
      }
    }
    return maxSurfDist;
  }
  
  
  /**
   * @brief RefineCurve
   * @param poly the curve to be refined
   * @param uniformFlag
   * 
   * Make one pass of refinement for all the edges of the curve that are distant from the basemesh
   * uses two parameters:
   * - par.minRefEdgeLen 
   * - par.surfDistThr
   */
    
  void RefineCurveByDistance(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    int startEdgeSize = poly.en;
    for(int i =0; i<startEdgeSize;++i)
    {
      EdgeType &ei = poly.edge[i];
      if(edge::Length(ei)>par.minRefEdgeLen)  
      {      
        CoordType farthestP, farthestN;
        ScalarType maxDist = MaxSegDist(ei.V(0),ei.V(1),farthestP, farthestN);
        if(maxDist > par.surfDistThr)  
        {
          edge::VEEdgeSplit(poly, &ei, farthestP, farthestN); 
        }
      }
    }
//    tri::Allocator<MeshType>::CompactEveryVector(poly);
//    printf("Refine %i -> %i\n",startEdgeSize,poly.en);fflush(stdout);
  }
  
  /**
   * @brief RefineCurveByBaseMesh
   * @param poly
   */
  
  void RefineCurveByBaseMesh(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    std::vector<int> edgeToRefineVec;
    for(int i=0; i<poly.en;++i) 
      edgeToRefineVec.push_back(i);
    int startEn=poly.en;  
    int iterCnt=0;
    while (!edgeToRefineVec.empty() && iterCnt<100) {
      iterCnt++;
      std::vector<int> edgeToRefineVecNext;
      for(int i=0; i<edgeToRefineVec.size();++i)
      {
        EdgeType &e = poly.edge[edgeToRefineVec[i]];
        CoordType splitPt;
        if(TestSplitSegWithMeshAdapt(e.V(0),e.V(1),splitPt))  
        {
          edge::VEEdgeSplit(poly, &e, splitPt); 
          edgeToRefineVecNext.push_back(edgeToRefineVec[i]);
          edgeToRefineVecNext.push_back(poly.en-1);
        } 
      }
      tri::Allocator<MeshType>::CompactEveryVector(poly);
      swap(edgeToRefineVecNext,edgeToRefineVec);
       printf("RefineCurveByBaseMesh %i en -> %i en\n",startEn,poly.en); fflush(stdout);    
    }
//
    SimplifyNullEdges(poly);
    SimplifyMidFace(poly);
    SimplifyMidEdge(poly);
    SnapPolyline(poly);    
    printf("RefineCurveByBaseMesh %i en -> %i en\n",startEn,poly.en); fflush(stdout);    
  }
  
  
  /**
   * @brief LaplacianFunctor basic Laplacian smoothing functor
   *
   * It computes the desired position for each vertex as the average of its
   * current position and the positions of its 1-ring neighbors. Used as the
   * position functor in the SmoothProject function.
   */
  struct LaplacianFunctor {
    std::vector<CoordType> operator()(const MeshType &poly) const {
      std::vector<CoordType> posVec(poly.vn, CoordType(0,0,0));
      std::vector<int>       cntVec(poly.vn, 0);
      for(int i=0; i<poly.en; ++i)
        for(int j=0; j<2; ++j) {
          int vi = tri::Index(poly, poly.edge[i].V0(j));
          posVec[vi] += poly.edge[i].V1(j)->P();
          cntVec[vi] += 1;
        }
      for(int i=0; i<poly.vn; ++i)
        posVec[i] = (poly.vert[i].P() + posVec[i]) / ScalarType(cntVec[i]+1);
      return posVec;
    }
  };

  /**
   * @brief QualityDistanceFieldFunctor quality field based smoothing functor
   *
   * @param poly the input curve mesh
   * @param com the CurveOnManifold class itself to quick access closest face 
   * @param scale a scaling factor to control the step size of the movement
   * along the quality, if 0 it will be automatically set to 1/2 of the CoM parameter `par.maxMoveDelta` (default is 0)
   * @param smoothBlend a blending factor to control the influence of the
   * smoothing (default is 0.5)
   *
   * It compute the new position using the quality field of the mesh assuming
   * that it is a distance field sampled per vertices and that we would like to move toward the zero of the distance field. 
   * We use gradient of the quality field for the direction and the sign of the quality for the versus of the direction. 
   * We move of a quantity proportional to the quality value at the vertex.
   *
   */

  struct QualityDistanceFieldFunctor
  {

    CoM<MeshType> &com;
    ScalarType scale;
    ScalarType smoothBlend = 0.5;
    QualityDistanceFieldFunctor(CoM<MeshType> &_com, ScalarType _scale=0, ScalarType _smoothBlend = 0.5) : com(_com), scale(_scale), smoothBlend(_smoothBlend) {};

    std::vector<CoordType> operator()(const MeshType &poly) const
    {
      // Step 1: Compute smoothed position using Laplacian smoothing
      std::vector<CoordType> smoothPosVec(poly.vn, CoordType(0, 0, 0));
      std::vector<int> cntVec(poly.vn, 0);
      for (int i = 0; i < poly.en; ++i)
        for (int j = 0; j < 2; ++j)
        {
          int vi = tri::Index(poly, poly.edge[i].V0(j));
          smoothPosVec[vi] += poly.edge[i].V1(j)->P();
          cntVec[vi] += 1;
        }
      for (int i = 0; i < poly.vn; ++i)
        smoothPosVec[i] = (poly.vert[i].P() + smoothPosVec[i]) / ScalarType(cntVec[i] + 1);

      // Step 2: Compute field-based position using the gradient of the quality field and moving toward zero
      std::vector<CoordType> fieldPosVec(poly.vn, CoordType(0, 0, 0));
      for (const VertexType &v : poly.vert)
      {
        CoordType ip; 
        FacePointer f = com.GetClosestFaceIP(v.P(),ip);
        if(!f) {
          printf("Fail to get closest face for vertex at position (%f, %f, %f)\n", v.P().X(), v.P().Y(), v.P().Z()); fflush(stdout);
        }
        else {
          ScalarType q = f->V(0)->Q() * ip[0] + f->V(1)->Q() * ip[1] + f->V(2)->Q() * ip[2];        
          CoordType fieldDir = GradientScalarField(*f, f->V(0)->Q(),f->V(1)->Q(),f->V(2)->Q());
          
          fieldPosVec[tri::Index(poly, v)] = v.P() + fieldDir * scale * q; // Move towards higher quality (lower distance)
        }
      }
      std::vector<CoordType> PosVec(poly.vn, CoordType(0, 0, 0));
      for (int i = 0; i < poly.vn; ++i)
        PosVec[i] = (smoothPosVec[i] * smoothBlend + fieldPosVec[i] * (1.0 - smoothBlend));

      return PosVec;
    }
  };

  /**
   * @brief MoveAndProject
   * @param poly
   * @param iterNum
   * @param moveWeight    [0..1] blend toward the desired position returned by the functor
   * @param projectWeight [0..1] blend toward the closest point on the surface
   * @param desiredPos    functor: std::vector<CoordType>(const MeshType&)
   *
   * Generic version of SmoothProject: the per-vertex desired position is
   * supplied by a functor instead of being hard-coded as a Laplacian average.
   */
  template<typename PositionFunctor>
  void MoveAndProject(MeshType &poly, int iterNum, ScalarType moveWeight, ScalarType projectWeight,
                      PositionFunctor &desiredPos)
  {
    tri::RequireCompactness(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    assert(poly.en>0 && base.fn>0);
    for(int k=0;k<iterNum;++k)
    {
      if(k==iterNum-1) projectWeight=1;
      std::vector<CoordType> desired = desiredPos(poly);

      for(int i=0; i<poly.vn; ++i)
        if(!poly.vert[i].IsS())
        {
          // Clamp the movement towards the desired position by maxMoveDist
          CoordType delta = desired[i] - poly.vert[i].P();
          ScalarType deltaLen = delta.Norm();
          if(deltaLen > par.maxMoveDelta) {
            delta *= (par.maxMoveDelta / deltaLen);
            desired[i] = poly.vert[i].P() + delta;
          } 

          CoordType newP = poly.vert[i].P()*(1.0-moveWeight) + desired[i]*moveWeight;
          
          CoordType closestP;
          FaceType *f = GetClosestFacePoint(newP, closestP);
          assert(f);
          poly.vert[i].P() = newP*(1.0-projectWeight) +closestP*projectWeight;
          poly.vert[i].N() = f->N();
        }
      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      RefineCurveByDistance(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      Simplify(poly);
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      int dupVertNum = Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(dupVertNum) {
        tri::Allocator<MeshType>::CompactEveryVector(poly);
        tri::UpdateTopology<MeshType>::VertexEdge(poly);
      }
    }
  }

  void SmoothProject(MeshType &poly, int iterNum, ScalarType smoothWeight, ScalarType projectWeight)
  {
    LaplacianFunctor lapFunct;
    MoveAndProject(poly, iterNum, smoothWeight, projectWeight, lapFunct);
  }


class EdgePointPred
{
public:
  std::map<std::pair<CoordType,CoordType>, VertexPointer> &edgeToPolyVertMap;
    
    EdgePointPred(std::map<std::pair<CoordType,CoordType>, VertexPointer> &_edgeToPolyVertMap):edgeToPolyVertMap(_edgeToPolyVertMap){};
    bool operator()(face::Pos<FaceType> ep) const
    {
      CoordType p0 = ep.V()->P();
      CoordType p1 = ep.VFlip()->P();
      if (p0>p1) std::swap(p0,p1);
      return edgeToPolyVertMap.find(std::make_pair(p0,p1)) != edgeToPolyVertMap.end();           
    }
};

struct EdgePointSplit
{
public:
  std::map<std::pair<CoordType,CoordType>, VertexPointer> &edgeToPolyVertMap;
  
  EdgePointSplit(std::map<std::pair<CoordType,CoordType>, VertexPointer> &_edgeToPolyVertMap):edgeToPolyVertMap(_edgeToPolyVertMap){};
    void operator()(VertexType &nv, face::Pos<FaceType> ep)
    {
      CoordType p0 = ep.V()->P();
      CoordType p1 = ep.VFlip()->P();
      if (p0>p1) std::swap(p0,p1);        
      VertexPointer vp=edgeToPolyVertMap[std::make_pair(p0,p1)];
      assert(vp);
      nv.P()=vp->P();
      return;
    }
    Color4b WedgeInterp(Color4b &c0, Color4b &c1)
    {
        Color4b cc;
        cc.lerp(c0,c1,0.5f);
        return Color4b::Red;
    }
    TexCoord2f WedgeInterp(TexCoord2f &t0, TexCoord2f &t1)
    {
        TexCoord2f tmp;
        assert(t0.n()== t1.n());
        tmp.n()=t0.n();
        tmp.t()=(t0.t()+t1.t())/2.0;
        return tmp;
    }
};

};

} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_CURVE_ON_SURF_H
