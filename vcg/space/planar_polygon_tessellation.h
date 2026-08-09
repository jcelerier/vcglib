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

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <vcg/space/point2.h>
#include <vcg/space/point3.h>

namespace vcg {

/** \addtogroup space */
/*@{*/
namespace planar_polygon_detail {

inline long double Orient2D(const Point2d &a, const Point2d &b, const Point2d &c)
{
	return (static_cast<long double>(b.X()) - static_cast<long double>(a.X())) *
	           (static_cast<long double>(c.Y()) - static_cast<long double>(a.Y())) -
	       (static_cast<long double>(b.Y()) - static_cast<long double>(a.Y())) *
	           (static_cast<long double>(c.X()) - static_cast<long double>(a.X()));
}

inline bool PointInTriangle(
	const Point2d &point,
	const Point2d &a,
	const Point2d &b,
	const Point2d &c,
	long double winding,
	long double epsilon)
{
	return winding * Orient2D(a, b, point) >= -epsilon
		&& winding * Orient2D(b, c, point) >= -epsilon
		&& winding * Orient2D(c, a, point) >= -epsilon;
}

} // namespace planar_polygon_detail

/**
 * Triangulate one finite, simple 2D polygon with ear clipping.
 *
 * The input must be one boundary loop, without holes, repeated consecutive
 * points, or self-intersections. Output indices are local to points, retain
 * its winding, and are appended only on success. Collinear vertices are supported;
 * no zero-area triangles are emitted. False means the polygon is degenerate,
 * invalid, or could not be triangulated reliably, and leaves output unchanged.
 */
template <class POINT_CONTAINER>
bool TessellatePlanarPolygon2(const POINT_CONTAINER &points, std::vector<int> &output)
{
	using namespace planar_polygon_detail;
	const size_t count = points.size();
	if (count < 3)
		return false;

	// Work in translated double coordinates. Translation avoids losing the small
	// differences that define the polygon when coordinates have a large offset.
	std::vector<Point2d> p(count);
	const double originX = double(points[0][0]);
	const double originY = double(points[0][1]);
	double scale = 0;
	for (size_t i = 0; i < count; ++i) {
		const double x = double(points[i][0]) - originX;
		const double y = double(points[i][1]) - originY;
		if (!std::isfinite(x) || !std::isfinite(y))
			return false;
		p[i] = Point2d(x, y);
		scale = std::max(scale, std::max(std::abs(x), std::abs(y)));
	}
	if (!(scale > 0))
		return false;

	const long double scale2 = static_cast<long double>(scale) * scale;
	const long double epsilon = scale2 * 64 * std::numeric_limits<double>::epsilon();
	const long double edgeEpsilon2 = scale2 * 1024
		* std::numeric_limits<double>::epsilon() * std::numeric_limits<double>::epsilon();

	long double signedDoubleArea = 0;
	for (size_t i = 0; i < count; ++i) {
		const Point2d &a = p[i];
		const Point2d &b = p[(i + 1) % count];
		const long double dx = static_cast<long double>(b.X()) - static_cast<long double>(a.X());
		const long double dy = static_cast<long double>(b.Y()) - static_cast<long double>(a.Y());
		if (dx * dx + dy * dy <= edgeEpsilon2)
			return false;
		signedDoubleArea += static_cast<long double>(a.X()) * static_cast<long double>(b.Y()) -
		                    static_cast<long double>(a.Y()) * static_cast<long double>(b.X());
	}
	if (std::abs(signedDoubleArea) <= epsilon * count)
		return false;
	const long double winding = signedDoubleArea > 0 ? 1 : -1;

	std::vector<int> previous(count), next(count);
	std::vector<unsigned char> active(count, 1);
	for (size_t i = 0; i < count; ++i) {
		previous[i] = int((i + count - 1) % count);
		next[i] = int((i + 1) % count);
	}

	std::vector<int> triangles;
	triangles.reserve(3 * (count - 2));
	int current = 0;
	int activeCount = int(count);
	while (activeCount > 2) {
		bool foundEar = false;
		for (int attempts = 0; attempts < activeCount; ++attempts) {
			const int before = previous[current];
			const int after = next[current];
			if (winding * Orient2D(p[before], p[current], p[after]) > epsilon) {
				// Convexity alone is insufficient: another active vertex inside
				// this triangle means the replacement diagonal is not an ear.
				bool containsVertex = false;
				for (size_t candidate = 0; candidate < count; ++candidate) {
					if (!active[candidate] || int(candidate) == before
						|| int(candidate) == current || int(candidate) == after)
						continue;
					if (PointInTriangle(
							p[candidate], p[before], p[current], p[after], winding, epsilon)) {
						containsVertex = true;
						break;
					}
				}
				if (!containsVertex) {
					triangles.push_back(before);
					triangles.push_back(current);
					triangles.push_back(after);
					next[before] = after;
					previous[after] = before;
					active[current] = 0;
					current = after;
					--activeCount;
					foundEar = true;
					break;
				}
			}
			current = next[current];
		}
		if (!foundEar)
			return false;
	}

	if (triangles.size() != 3 * (count - 2))
		return false;

	// Validate the postcondition independently from ear selection. Consistent
	// winding and equal area reject overlaps, outside triangles, and gaps.
	long double triangleDoubleArea = 0;
	for (size_t i = 0; i < triangles.size(); i += 3) {
		const long double area = winding * Orient2D(
			p[triangles[i]], p[triangles[i + 1]], p[triangles[i + 2]]);
		if (area <= epsilon)
			return false;
		triangleDoubleArea += area;
	}
	const long double areaTolerance = std::max(
		std::abs(signedDoubleArea) * 1e-12L, epsilon * count * 4);
	if (std::abs(triangleDoubleArea - std::abs(signedDoubleArea)) > areaTolerance)
		return false;

	output.insert(output.end(), triangles.begin(), triangles.end());
	return true;
}

/**
 * Triangulate one finite, simple, approximately planar polygon embedded in 3D.
 *
 * A Newell area vector supplies a deterministic normal. Coordinates relative
 * to the first point are projected to the dominant coordinate plane, avoiding
 * an unstable arbitrary basis and large-offset cancellation. The 2D routine
 * defines output and failure behavior. False also reports a degenerate normal
 * or a departure from planarity larger than the input scalar precision allows.
 */
template <class POINT_CONTAINER>
bool TessellatePlanarPolygon3(const POINT_CONTAINER &points, std::vector<int> &output)
{
	typedef typename POINT_CONTAINER::value_type Point3x;
	typedef typename Point3x::ScalarType ScalarType;
	const size_t count = points.size();
	if (count < 3)
		return false;

	const Point3d origin{
		static_cast<double>(points[0][0]),
		static_cast<double>(points[0][1]),
		static_cast<double>(points[0][2])};
	std::vector<Point3d> relative(count);
	double scale = 0;
	for (size_t i = 0; i < count; ++i) {
		relative[i] = Point3d(
			double(points[i][0]) - origin.X(),
			double(points[i][1]) - origin.Y(),
			double(points[i][2]) - origin.Z());
		if (!std::isfinite(relative[i].X()) || !std::isfinite(relative[i].Y())
			|| !std::isfinite(relative[i].Z()))
			return false;
		scale = std::max(scale, relative[i].Norm());
	}
	if (!(scale > 0))
		return false;

	Point3d normal(0, 0, 0);
	for (size_t i = 0; i < count; ++i)
		normal += relative[i] ^ relative[(i + 1) % count];
	const double normalNorm = normal.Norm();
	const double normalTolerance = scale * scale * count * 64
		* std::numeric_limits<double>::epsilon();
	if (!(normalNorm > normalTolerance))
		return false;
	normal /= normalNorm;

	double scalarEpsilon = double(std::numeric_limits<ScalarType>::epsilon());
	if (!(scalarEpsilon > 0))
		scalarEpsilon = std::numeric_limits<double>::epsilon();
	const double planarityTolerance = scale * std::max(1e-12, 32 * scalarEpsilon);
	for (size_t i = 0; i < count; ++i)
		if (std::abs(relative[i] * normal) > planarityTolerance)
			return false;

	std::vector<Point2d> projected;
	projected.reserve(count);
	const Point3d absoluteNormal(std::abs(normal.X()), std::abs(normal.Y()), std::abs(normal.Z()));
	for (size_t i = 0; i < count; ++i) {
		if (absoluteNormal.X() >= absoluteNormal.Y() && absoluteNormal.X() >= absoluteNormal.Z())
			projected.push_back(Point2d(relative[i].Y(), relative[i].Z()));
		else if (absoluteNormal.Y() >= absoluteNormal.Z())
			projected.push_back(Point2d(relative[i].X(), relative[i].Z()));
		else
			projected.push_back(Point2d(relative[i].X(), relative[i].Y()));
	}
	return TessellatePlanarPolygon2(projected, output);
}

/*@}*/
} // end namespace
#endif
