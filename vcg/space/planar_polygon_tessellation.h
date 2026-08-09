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

#ifndef __VCGLIB_PLANAR_POLYGON_TESSELLATOR
#define __VCGLIB_PLANAR_POLYGON_TESSELLATOR

#include <assert.h>
#include <vector>
#include <vcg/space/segment2.h>
#include <vcg/space/point3.h>
#include <vcg/math/random_generator.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
    /**
     * Given four 2D points p00,p01,p10,p11, it returns true if the segments p00-p01 and p10-p11 intersect.
	*/
	template <class ScalarType> 
	bool Cross(	 const Point2<ScalarType> & p00,
				 const Point2<ScalarType> & p01,
				 const Point2<ScalarType> & p10,
				 const Point2<ScalarType> & p11)
	{
		Point2<ScalarType> vec0 = p01-p00;
		Point2<ScalarType> vec1 = p11-p10;
		if ( ( vec0^ (p11-p00)) *  ( vec0^ (p10 - p00)) >=0) return false;
		if ( ( vec1^ (p01-p10)) *  ( vec1^ (p00 - p10)) >=0) return false;
		return true;
	}

	template <class S>
	bool Intersect(size_t cur , int v2, std::vector<int> & next, std::vector<Point2<S> > & points2){
		for(size_t i  = 0; i < points2.size();++i)
			if( (next[i]!=-1) && (i!=cur))
				if( Cross(points2[cur], points2[v2],points2[i],points2[next[i]]))
					return true;
		return false;
	}


	/**
	 * Triangulate one simple 2D polygon with ear clipping.
	 *
	 * The input points describe a single boundary loop in either orientation;
	 * holes and self-intersections are not supported. Output indices refer to
	 * positions in points2 and triangles retain the input winding. An ear is
	 * accepted when it is convex and its replacement diagonal does not cross
	 * any active boundary edge. Degenerate input falls back to a fan over the
	 * remaining loop so callers still receive a complete n-2 triangulation.
	 *
	 * output is appended to rather than cleared.
	 */
	template <class POINT_CONTAINER>
	void TessellatePlanarPolygon2( POINT_CONTAINER &  points2, std::vector<int> & output){
		typedef typename POINT_CONTAINER::value_type Point2x;
		typedef typename Point2x::ScalarType S;
		if(points2.size() < 3) return;
		// tessellate
		//  first very inefficient implementation
		std::vector<int> next,prev;
		for(size_t i = 0; i < points2.size(); ++i) next.push_back((i+1)%points2.size());
		for(size_t i = 0; i < points2.size(); ++i) prev.push_back((i+points2.size()-1)%points2.size());
		int v1,v2;
		// check orientation
		S orient = 0.0;
		for(size_t i = 0 ; i < points2.size(); ++i)
			orient += points2[i] ^ points2[next[i]];
		orient = (orient>=0)? 1.0:-1.0;

		int cur = 0;
		int activeCnt = int(points2.size());
		while(activeCnt > 2)
		{
			bool earFound = false;
			int attempts = 0;
			while(attempts < activeCnt)
			{
				v1 = next[cur];
				v2 = next[v1];
				// Removing v1 replaces boundary edges cur-v1 and v1-v2 with
				// cur-v2. Convexity plus a non-crossing replacement diagonal
				// identifies a valid ear for a simple polygon.
				if( ( (orient*((points2[v1] - points2[cur]) ^ (points2[v2] - points2[cur]))) >= 0.0) &&
					!Intersect(cur, v2,next,points2))
				{
					// output the face
					output.push_back(cur);
					output.push_back(v1);
					output.push_back(v2);

					// readjust the topology: remove v1 from active ring
					next[cur] = v2;
					prev[v2] = cur;
					prev[v1] = -1;
					next[v1] = -1;
					--activeCnt;
					earFound = true;
					break;
				}
				do{cur = (cur+1)%points2.size();} while(next[cur]==-1);
				++attempts;
			}

			if(earFound)
			{
				do{cur = (cur+1)%points2.size();} while(next[cur]==-1);
				continue;
			}

			// Fallback for degenerate / numerically problematic cases:
			// triangulate remaining active ring as a fan to guarantee n-2 output triangles.
			std::vector<int> ring;
			ring.reserve(size_t(activeCnt));
			int start = -1;
			for(size_t i = 0; i < next.size(); ++i)
				if(next[i] != -1) { start = int(i); break; }
			if(start == -1)
				break;
			int it = start;
			do {
				ring.push_back(it);
				it = next[it];
			} while(it != -1 && it != start && ring.size() <= size_t(activeCnt));

			if(ring.size() < 3)
				break;
			for(size_t i = 1; i + 1 < ring.size(); ++i)
			{
				output.push_back(ring[0]);
				output.push_back(ring[i]);
				output.push_back(ring[i+1]);
			}
			break;
		}
	}

	/**
	 * Triangulate one simple planar polygon embedded in 3D.
	 *
	 * A stable projection plane is selected from a large-area input triangle,
	 * the polygon is projected to 2D, and TessellatePlanarPolygon2 performs the
	 * ear clipping. Input vertices must form one approximately planar boundary;
	 * holes and self-intersections are not supported. Output is a flat sequence
	 * of triangle indices local to points, preserving the boundary winding, and
	 * is appended to rather than cleared.
	 */

	template <class POINT_CONTAINER>
	void TessellatePlanarPolygon3( POINT_CONTAINER &  points, std::vector<int> & output){
		typedef typename POINT_CONTAINER::value_type Point3x;
		typedef typename Point3x::ScalarType S;
		Point3x n;
		if(points.size()==3)
		{
			output.push_back(0);
			output.push_back(1);
			output.push_back(2);
			return;
		}
		
		// if the polygon is a quad we can optimize the tessellation just checking to the shortest diagonal
		if(points.size()==4)
		{
			if(Distance(points[0],points[2])<Distance(points[1],points[3])){
				output.push_back(0);
				output.push_back(1);
				output.push_back(2);
				output.push_back(0);
				output.push_back(2);
				output.push_back(3);
			}else{
				output.push_back(0);
				output.push_back(1);
				output.push_back(3);
				output.push_back(1);
				output.push_back(2);
				output.push_back(3);
			}
			return;
		}
		
		math::SubtractiveRingRNG rg;
		size_t i12[2]={0,0};
		S bestsn = -1.0;
		Point3x bestn,u,v;
		
		// find the best normal for projection on the plane
		for(size_t i  =0; i < points.size();++i){
			do{
				i12[0] = rg.generate(points.size());
			} while(i12[0]==i);
			do{
				i12[1] = rg.generate(points.size());
			} while((i12[1]==i || i12[1]==i12[0]));
			if(!(i12[0]!=i12[1] && i12[0]!=i))
				assert(0);
			n = (points[i12[0]]-points[i])^(points[i12[1]]-points[i]);
			S sn = n.SquaredNorm();
			if(sn > bestsn){ bestsn = sn; bestn = n;} 
		}
		
		GetUV(bestn,u,v);
		// project the coordinates
		std::vector<Point2<S> > points2;
		for(size_t i = 0; i < points.size(); ++i){
			Point3x & p = points[i];
			points2.push_back(Point2<S>(p*u,p*v));
		}
		TessellatePlanarPolygon2( points2,output);
	}

/*@}*/
} // end namespace
#endif
