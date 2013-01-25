/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_COLLIDE_2D_SPHERE_CONVEX_H
#define BT_COLLIDE_2D_SPHERE_CONVEX_H

inline btScalar btPointLineDistance3d(const btVector3& lineA, const btVector3& lineB, const btVector3& point)
{
	//Line equation: lineA + lineAToB * t
	btVector3 lineAToB = (lineB - lineA).normalized();
	
	btVector3 pointToLineA = lineA - point;

	btVector3 pointToLineNearest = pointToLineA - pointToLineA.dot(lineAToB)*lineAToB;
	
	return pointToLineNearest.length();
}

//Input: 
//	Point
//	2D line segment
//Output:
//	Normal on line segment
//	Nearest point on line segment
//	Distance from point to line segment
//	Normal on line
inline void btPointLineSegmentCollision2d(btScalar lineA_u, btScalar lineA_v, btScalar lineB_u, btScalar lineB_v, 
										btScalar point_u, btScalar point_v, 
										btScalar& out_normalOnSegment_u, btScalar& out_normalOnSegment_v,
										btScalar& out_pointOnSegment_u, btScalar& out_pointOnSegment_v, btScalar& out_segmentDistance,
										btScalar& out_normalOnLine_u, btScalar& out_normalOnLine_v)
{
	btScalar lineAToB_u = lineB_u - lineA_u;
	btScalar lineAToB_v = lineB_v - lineA_v;
	
	btScalar pointToLineA_u = lineA_u - point_u;
	btScalar pointToLineA_v = lineA_v - point_v;
	
	btScalar lineAToB_length2 = lineAToB_u*lineAToB_u + lineAToB_v*lineAToB_v;
	
	btScalar t_line = -(pointToLineA_u * lineAToB_u + pointToLineA_v * lineAToB_v) / lineAToB_length2;
	
	//Line segment: return normal, nearest point, and distance
	{
		btScalar t_lineSegment = btMin( btMax(btScalar(0.0), t_line), btScalar(1.0) );
		
		btScalar nearestPointOnSegment_u = lineA_u + lineAToB_u * t_lineSegment; 
		btScalar nearestPointOnSegment_v = lineA_v + lineAToB_v * t_lineSegment; 
		
		btScalar segmentToPoint_u = point_u - nearestPointOnSegment_u;
		btScalar segmentToPoint_v = point_v - nearestPointOnSegment_v;
		
		btScalar segmentToPoint_length = btSqrt(segmentToPoint_u*segmentToPoint_u + segmentToPoint_v*segmentToPoint_v);
		
		out_normalOnSegment_u = (segmentToPoint_length > SIMD_EPSILON) ? segmentToPoint_u / segmentToPoint_length : btScalar(0.0);
		out_normalOnSegment_v = (segmentToPoint_length > SIMD_EPSILON) ? segmentToPoint_v / segmentToPoint_length : btScalar(1.0);
		
		out_segmentDistance = segmentToPoint_length;
		out_pointOnSegment_u = nearestPointOnSegment_u;
		out_pointOnSegment_v = nearestPointOnSegment_v;
	}
	
	//Line: return normal on line
	{
		btScalar nearestPointOnLine_u = lineA_u + lineAToB_u * t_line; 
		btScalar nearestPointOnLine_v = lineA_v + lineAToB_v * t_line; 
		
		btScalar lineToPoint_u = point_u - nearestPointOnLine_u;
		btScalar lineToPoint_v = point_v - nearestPointOnLine_v;
		
		btScalar lineToPoint_length = btSqrt(lineToPoint_u*lineToPoint_u + lineToPoint_v*lineToPoint_v);
		
		out_normalOnLine_u = (lineToPoint_length > SIMD_EPSILON) ? lineToPoint_u / lineToPoint_length : btScalar(0.0);
		out_normalOnLine_v = (lineToPoint_length > SIMD_EPSILON) ? lineToPoint_v / lineToPoint_length : btScalar(1.0);
	}
}

///Detects collisions between a sphere(circle) and convex polygon in 2D.
///polygonEdges_u[], and polygonEdges_u[] should be arrays with numEdges+1 elements, and should enclose polygonCenter_u/v.
///Edges are expressed as a series of connected line segments; for instance, polygonEdges_u[0] and polygonEdges_u[1]
///are the u coordinates of the first line segment, polygonEdges_u[1] and polygonEdges_u[2] are the coordinates of the second, and so on.
inline void btCollide2dSphereConvex(btScalar spherePos_u, btScalar spherePos_v, btScalar sphereRadius, 
							btScalar polygonCenter_u, btScalar polygonCenter_v, btScalar* polygonEdges_u, btScalar* polygonEdges_v, int numEdges, 
							bool& out_sphereInPolygon, btScalar& out_distance2d,
							btScalar& out_contact2d_u, btScalar& out_contact2d_v,
							btScalar& out_normal2d_u, btScalar& out_normal2d_v)
{	
	btScalar nearestLineSegmentDistance(BT_LARGE_FLOAT);
	btScalar nearestLineSegmentContact_u, nearestLineSegmentContact_v;
	btScalar nearestLineSegmentNormal_u, nearestLineSegmentNormal_v;
			
	bool sphereInPolygon = true;
	for(int i = 0; i < numEdges; ++i)
	{
		btScalar lineA_u = polygonEdges_u[i];
		btScalar lineA_v = polygonEdges_v[i];
		
		btScalar lineB_u = polygonEdges_u[i+1];
		btScalar lineB_v = polygonEdges_v[i+1];
		
		//Calculate normal by projecting the polygon center onto the line
		btScalar discard[5];	//Unused
		btScalar rectangleNormal_u, rectangleNormal_v;
		btPointLineSegmentCollision2d(lineA_u, lineA_v, lineB_u, lineB_v, polygonCenter_u, polygonCenter_v, 
									discard[0], discard[1], discard[2], discard[3], discard[4],
									rectangleNormal_u, rectangleNormal_v);
		rectangleNormal_u = -rectangleNormal_u;	//Since the polygon center is used, the normal will point inwards
		rectangleNormal_v = -rectangleNormal_v;	//Invert it to make it point outwards
		
		//Use normals to test which side of line sphere is on
		btScalar normalOnSegment_u, normalOnSegment_v, pointOnSegment_u, pointOnSegment_v, distanceToSegment; 
		btScalar normalOnLine_u, normalOnLine_v; 
		btPointLineSegmentCollision2d(lineA_u, lineA_v, lineB_u, lineB_v, spherePos_u, spherePos_v, 
									normalOnSegment_u, normalOnSegment_v, pointOnSegment_u, pointOnSegment_v, distanceToSegment,
									normalOnLine_u, normalOnLine_v);
		
		btScalar dotProduct = rectangleNormal_u*normalOnLine_u + rectangleNormal_v*normalOnLine_v;
		
		
		if(dotProduct > 0.0) sphereInPolygon = false; 	//Sphere is not on the same side of the line as the polygon center
		
		if(distanceToSegment < nearestLineSegmentDistance)
		{
			nearestLineSegmentDistance = distanceToSegment;
			nearestLineSegmentContact_u = pointOnSegment_u;
			nearestLineSegmentContact_v = pointOnSegment_v;
			nearestLineSegmentNormal_u = normalOnSegment_u;
			nearestLineSegmentNormal_v = normalOnSegment_v;
		}
	}
	
	//
	out_sphereInPolygon = sphereInPolygon;
	out_contact2d_u = nearestLineSegmentContact_u;
	out_contact2d_v = nearestLineSegmentContact_v;
	
	if(!sphereInPolygon)
	{
		out_distance2d = nearestLineSegmentDistance - sphereRadius;
		
		//Since the sphere is outside the polygon, the line normal already points outwards
		out_normal2d_u = nearestLineSegmentNormal_u;
		out_normal2d_v = nearestLineSegmentNormal_v;
	}
	else
	{
		out_distance2d = -nearestLineSegmentDistance - sphereRadius;
	
		//The sphere is inside the polygon, so the line normal will point into the polygon
		//Invert the line segment normal to make it point away from the polygon center
		out_normal2d_u = -nearestLineSegmentNormal_u;
		out_normal2d_v = -nearestLineSegmentNormal_v;
	}
}

#endif