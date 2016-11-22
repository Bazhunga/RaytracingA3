/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur. This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.


	// The ray is shot out in camera coordinates. But remember that
	// we had previously changed that to world coordinates! 
	// We need to change this back to camera coordinates to make 
	// things easier 
	// We need to take this and turn it into world coordinates by
	// multiplying it with modeltoWorld
	Ray3D object_space_ray; 
	object_space_ray.origin = worldToModel * ray.origin;
	object_space_ray.dir = worldToModel * ray.dir;

	// Now we're in camera coordinates
	// We want to see if it intersects with the unit square
	// Since this lies on the xy plane, we just see when it hits z = 0

	// Ray is p(t) = q + rt where q is the ray.origin and r is ray.dir (unit vector)
	// Just look at the Z element.
	// q_z + r_z*t = 0 
	// If t is negative, there is no intersection (we're shooting off to the other direction)
	// If t is positive, then do a check
	// Check if the x and y coordinates lie between the unit square


	double ray_test_t = -object_space_ray.origin[2]/object_space_ray.dir[2];

	if (ray_test_t > 0) {
		double plane_x = object_space_ray.origin[0]+object_space_ray.dir[0]*ray_test_t;
		double plane_y = object_space_ray.origin[1]+object_space_ray.dir[1]*ray_test_t;

		if (plane_x<=0.5 && plane_x>=-0.5 && plane_y<=0.5 && plane_y>=-0.5) {
			Point3D intersection_point(plane_x, plane_y, 0);
			ray.intersection.point = modelToWorld * intersection_point;
			Vector3D normal(0, 0, 1);
			ray.intersection.normal = modelToWorld * normal;
			ray.intersection.t_value = ray_test_t;
			ray.intersection.none = false;

			return true;
		}
	}
	ray.intersection.none = true;
	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	// Ray3D object_space_ray 
	// object_space_ray.origin = worldToModel * ray.origin;
	// object_space_ray.dir = worldToModel * ray.dir;

	// Intersect with the unit sphere centered around the origin 
	// Draw a vector from center of circle (0,0,0) to the ray origin 
	// Project that onto yo


	
	return false;
}

 