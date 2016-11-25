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
	Ray3D object_space_ray = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);

	// Now we're in camera coordinates
	// We want to see if it intersects with the unit square
	// Since this lies on the xy plane, we just see when it hits z = 0

	// Ray is p(t) = q + rt where q is the ray.origin and r is ray.dir (unit vector)
	// Just look at the Z element.
	// q_z + r_z*t = 0 
	// If t is negative, there is no intersection (we're shooting off to the other direction)
	// If t is positive, then do a check
	// Check if the x and y coordinates lie between the unit square


	double ray_test_t = -(object_space_ray.origin[2]/object_space_ray.dir[2]);

	Point3D plane_intersect_point = object_space_ray.origin + ray_test_t * object_space_ray.dir;
	// std::cout << "plane_intersect_point: "<< plane_intersect_point[0]<< " " << plane_intersect_point[1] << " " << plane_intersect_point[2] << ".\n";
	if (ray_test_t > 0) {
		if (plane_intersect_point[0]<=0.5 && plane_intersect_point[0]>=-0.5 
			&& plane_intersect_point[1]<=0.5 && plane_intersect_point[1]>=-0.5) {
			if(ray.intersection.none || ray_test_t < ray.intersection.t_value){
				Intersection intersection_point;
				intersection_point.point = modelToWorld * plane_intersect_point;
				Vector3D normal(0, 0, 1);
				intersection_point.normal = worldToModel.transpose() * normal;
				intersection_point.t_value = ray_test_t;
				intersection_point.none = false;
				ray.intersection = intersection_point;
				return true;
			}
		}
	}
	// Do not reset ray.intersection.non to true!! It might have gone through the 
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
	Ray3D object_space_ray = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);

	object_space_ray.dir.normalize(); // We need to normalize to get the correct 

   // q - c where c = (0, 0, 0)
   Vector3D ray_to_circle_center(-1*object_space_ray.origin[0], -1*object_space_ray.origin[1], -1*object_space_ray.origin[2]);

   // y^2 = |q - c|^2 - |(q - c) . r|^2
   double ray_test_y_squared = pow(ray_to_circle_center.length(),2) - pow(ray_to_circle_center.dot(object_space_ray.dir),2);

   // k^2 = d^2 - y^2 where d = 1
   double ray_test_k_squared = 1 - ray_test_y_squared;
   // if (ray_test_k_squared > 0) {
   // 	std::cout << "K found: "<< ray_test_k_squared << ".\n";
   // }

   // if k exists, there is at least one point of intersection
	if (ray_test_k_squared >= 0) {
		double k = sqrt(ray_test_k_squared);

		Vector3D ray_including_k;
		Vector3D ray_to_intersection;

		ray_including_k = object_space_ray.dir.dot(ray_to_circle_center) * object_space_ray.dir;

		// multiply with ray_including_k to reduce length to length - k
		double length_factor = (ray_including_k.length() - k);
		ray_to_intersection = length_factor * object_space_ray.dir;

	    double t_value = (ray_to_intersection.length() / object_space_ray.dir.length());

	   	if (ray.intersection.none || t_value < ray.intersection.t_value){

	      // (x,y,z) = r + q
	      double x = ray_to_intersection[0] + object_space_ray.origin[0];
	      double y = ray_to_intersection[1] + object_space_ray.origin[1];
	      double z = ray_to_intersection[2] + object_space_ray.origin[2];

	      Point3D intersection_point(x, y, z);

	      ray.intersection.point = modelToWorld * intersection_point;


	      // Normal vector at the point of intersection
	      Vector3D normal(intersection_point[0], intersection_point[1], intersection_point[2]);
	      ray.intersection.normal = worldToModel.transpose() * normal;
	      ray.intersection.normal.normalize();

	      // All the t_value variables should have the same value
	      ray.intersection.t_value = t_value;

	      ray.intersection.none = false;
	      return true;

	   	}
	}
	return false;
}

 