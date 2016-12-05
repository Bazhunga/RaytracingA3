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
	object_space_ray.intersection = ray.intersection;
	std::string previousIntersection = ray.previousShape;

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
			// Check if previous intersection is square. If yes, ignore this. 
			if (previousIntersection.compare("square") != 0){ 
				if((ray.intersection.none || ray_test_t < ray.intersection.t_value)){
					Intersection intersection_point;
					intersection_point.point = modelToWorld * plane_intersect_point;
					Vector3D normal(0, 0, 1);
					if(object_space_ray.dir.dot(normal) > 1) {
						normal = -1 * normal; // Normal was facing wrong way
					}
					intersection_point.normal = worldToModel.transpose() * normal;
					intersection_point.t_value = ray_test_t;
					intersection_point.none = false;
					ray.intersection = intersection_point;
					ray.intersection.normal.normalize();
					ray.intersection.shape = "square";
					return true;
				}
			}
			else {
				// We're self intersecting
				// Reset the ray.intersection.none because this is a reflected ray
				ray.intersection.none = true; 
			}
		}
	}
	// Do not reset ray.intersection.non to true!! It might have gone through the circle 
	return false;
}

bool UnitSquare2::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
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
	object_space_ray.intersection = ray.intersection;
	std::string previousIntersection = ray.previousShape;

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
			// Check if previous intersection is square. If yes, ignore this. 
			if (previousIntersection.compare("square2") != 0){ 
				if((ray.intersection.none || ray_test_t < ray.intersection.t_value)){
					Intersection intersection_point;
					intersection_point.point = modelToWorld * plane_intersect_point;
					Vector3D normal(0, 0, 1);
					if(object_space_ray.dir.dot(normal) > 1) {
						normal = -1 * normal; // Normal was facing wrong way
					}
					intersection_point.normal = worldToModel.transpose() * normal;
					intersection_point.t_value = ray_test_t;
					intersection_point.none = false;
					ray.intersection = intersection_point;
					ray.intersection.normal.normalize();
					ray.intersection.shape = "square2";
					return true;
				}
			}
			else {
				// We're self intersecting
				// Reset the ray.intersection.none because this is a reflected ray
				ray.intersection.none = true; 
			}
		}
	}
	// Do not reset ray.intersection.non to true!! It might have gone through the circle 
	return false;
}

bool UnitSquareTextured::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
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
	object_space_ray.intersection = ray.intersection;
	std::string previousIntersection = ray.previousShape;


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
			if (previousIntersection.compare("square") != 0){ 
				if((ray.intersection.none || ray_test_t < ray.intersection.t_value)){
					Intersection intersection_point;
					intersection_point.point = modelToWorld * plane_intersect_point;
					Vector3D normal(0, 0, 1);
					intersection_point.normal = worldToModel.transpose() * normal;
					intersection_point.t_value = ray_test_t;
					intersection_point.none = false;
					ray.intersection = intersection_point;
					ray.intersection.normal.normalize();
					ray.intersection.shape = "square";
					return true;
				}
			}
			else {
				// We're self intersecting
				// Reset the ray.intersection.none because this is a reflected ray
				ray.intersection.none = true; 
			}
		}
	}
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
	Vector3D osr_dir_unit = object_space_ray.dir;
	osr_dir_unit.normalize();

	object_space_ray.intersection = ray.intersection;
	// object_space_ray.dir.normalize(); // We need to normalize to get the correct 
	std::string previousIntersection = ray.previousShape;

   // q - c where c = (0, 0, 0)
   Vector3D ray_to_circle_center(-1*object_space_ray.origin[0], -1*object_space_ray.origin[1], -1*object_space_ray.origin[2]);
   double rtc = ray_to_circle_center.length();
   double proj = ray_to_circle_center.dot(osr_dir_unit);

   // y^2 = |q - c|^2 - |(q - c) . r|^2
   double ray_test_y_squared = pow(rtc,2) - pow(proj,2);

   // k^2 = d^2 - y^2 where d = 1
   double ray_test_k_squared = 1 - ray_test_y_squared;
   // if (ray_test_k_squared > 0) {
   // 	std::cout << "K found: "<< ray_test_k_squared << ".\n";
   // }

   // if k exists, there is at least one point of intersection
	if (ray_test_k_squared >= 0) {
		double k = sqrt(ray_test_k_squared);

		Vector3D ray_to_intersection;

		// multiply with ray_including_k to reduce length to length - k
		ray_to_intersection = (proj - k) * osr_dir_unit;

	   double t_value = (ray_to_intersection.length() / object_space_ray.dir.length());
	   if(previousIntersection.compare("sphere") != 0) {
	   	if (ray.intersection.none || t_value < ray.intersection.t_value){

		      // (x,y,z) = r + q
		      double x = ray_to_intersection[0] + object_space_ray.origin[0];
		      double y = ray_to_intersection[1] + object_space_ray.origin[1];
		      double z = ray_to_intersection[2] + object_space_ray.origin[2];

		      Point3D intersection_point(x, y, z);

		      ray.intersection.point = modelToWorld * intersection_point;


		      // Normal vector at the point of intersection
		      Vector3D normal(intersection_point[0], intersection_point[1], intersection_point[2]);
		      if(object_space_ray.dir.dot(normal) > 1) {
					normal = -1 * normal; // Normal was facing wrong way
				}
		      // ray.intersection.normal = worldToModel.transpose() * normal;
		      ray.intersection.normal = worldToModel.transpose() * normal;
		      ray.intersection.normal.normalize();
		      // ray.intersection.normal.normalize();

		      // All the t_value variables should have the same value
		      ray.intersection.t_value = t_value;

		      ray.intersection.none = false;
		     	ray.intersection.shape = "sphere";
		      return true;
	   	}
	   }
	   else {
	      // ray.intersection.none = true;
	   }
	}
	return false;
}

bool UnitSphere2::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	Ray3D object_space_ray = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);
	Vector3D osr_dir_unit = object_space_ray.dir;
	osr_dir_unit.normalize();

	object_space_ray.intersection = ray.intersection;
	std::string previousIntersection = ray.previousShape;

   Vector3D ray_to_circle_center(-1*object_space_ray.origin[0], -1*object_space_ray.origin[1], -1*object_space_ray.origin[2]);
   double rtc = ray_to_circle_center.length();
   double proj = ray_to_circle_center.dot(osr_dir_unit);

   double ray_test_y_squared = pow(rtc,2) - pow(proj,2);
   double ray_test_k_squared = 1 - ray_test_y_squared;
	if (ray_test_k_squared >= 0) {
		double k = sqrt(ray_test_k_squared);

		Vector3D ray_to_intersection;

		ray_to_intersection = (proj - k) * osr_dir_unit;

	   double t_value = (ray_to_intersection.length() / object_space_ray.dir.length());
	   if(previousIntersection.compare("sphere2") != 0) {
	   	if (ray.intersection.none || t_value < ray.intersection.t_value){

		      double x = ray_to_intersection[0] + object_space_ray.origin[0];
		      double y = ray_to_intersection[1] + object_space_ray.origin[1];
		      double z = ray_to_intersection[2] + object_space_ray.origin[2];

		      Point3D intersection_point(x, y, z);

		      ray.intersection.point = modelToWorld * intersection_point;

		      Vector3D normal(intersection_point[0], intersection_point[1], intersection_point[2]);
		      if(object_space_ray.dir.dot(normal) > 1) {
					normal = -1 * normal; // Normal was facing wrong way
				}
		      // ray.intersection.normal = worldToModel.transpose() * normal;
		      ray.intersection.normal = worldToModel.transpose() * normal;
		      ray.intersection.normal.normalize();
		      // ray.intersection.normal.normalize();

		      // All the t_value variables should have the same value
		      ray.intersection.t_value = t_value;

		      ray.intersection.none = false;
		     	ray.intersection.shape = "sphere2";
		      return true;
	   	}
	   }
	   else {
	      ray.intersection.none = true;
	   }
	}
	return false;
}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// // Unit cylinder: circle of radius 1 along the y axis between 0.5 and -0.5

	// // Transform ray into object space
	// Point3D origin = worldToModel * ray.origin;
	// Vector3D dir = worldToModel * ray.dir;
	// Ray3D object_space_ray = Ray3D(origin, dir);

	// double a = pow(object_space_ray.dir[0],2) + pow(object_space_ray.dir[2],2);
	// double b = 2*object_space_ray.dir[0]*object_space_ray.origin[0] + 2 * object_space_ray.dir[2] * object_space_ray.origin[2];
	// double c = pow(object_space_ray.origin[0],2) + pow(object_space_ray.origin[2],2) - 1;

	// double discrim = b*b - 4*a*c;

	// if (discrim >= 0) {
	// 	double t0 = ((-1)*b + sqrt(discrim))/(2*a);
	// 	double t1 = ((-1)*b - sqrt(discrim))/(2*a);

	// 	double t = -1;

	// 	// Reordering
	// 	if (t0 > t1) {
	// 		double tmp = t0;
	// 		t0 = t1;
	// 		t1 = tmp;
	// 	}

	// 	// Get the y values
	// 	double y0 = origin[1] + t0 * dir[0];
	// 	double y1 = origin[1] + t1 * dir[0];

	// 	if ((y0 > 1 && y1 > 1) || (y0 < -1 && y1 < -1)){ // total wiff
	// 		return false;
	// 	}

	// 	if (y0 > 1 && y1 < 1) { // Hit cylinder cap at +1			
	// 		double t_cap = t0 + (t1-t0) * (y0-1) / (y0-y1);
	// 		if (t_cap<=0) {
	// 			return false;
	// 		} 
	// 		if (ray.intersection.none || t_cap < ray.intersection.t_value) {
	// 			std::cout<< "yas1\n";
	// 			Point3D intersection_point = origin + t_cap * dir;
	// 			ray.intersection.t_value = t_cap;
	// 			ray.intersection.point =  modelToWorld * intersection_point;
	// 			std::cout << "X: " << ray.intersection.point[0] << " Y: " << ray.intersection.point[1] << " Z: "<< ray.intersection.point[2] << "\n";
	// 			ray.intersection.normal = worldToModel.transpose() * Vector3D(0, 1, 0);
	// 			ray.intersection.normal.normalize();
	// 			ray.intersection.none = false;
	// 			ray.intersection.shape = "cylinder";
	// 			return true;
	// 		}
	// 	}

	// 	if (y0 < -1 && y1 > -1) { // Hit cylinder cap at -1
	// 		double t_cap = t0 + (t1-t0) * (y0-1) / (y0-y1);
	// 		if (t_cap<=0) {
	// 			return false;
	// 		} 
	// 		if (ray.intersection.none || t_cap < ray.intersection.t_value) {
	// 			std::cout<< "yas2\n";
	// 			Point3D intersection_point = origin + t_cap * dir;
	// 			ray.intersection.t_value = t_cap;
	// 			ray.intersection.point =  modelToWorld * intersection_point;
	// 			ray.intersection.normal = worldToModel.transpose() * Vector3D(0, -1, 0);
	// 			ray.intersection.normal.normalize();
	// 			ray.intersection.none = false;
	// 			ray.intersection.shape = "cylinder";
	// 			return true;
	// 		}
	// 	}

	// 	if (y0 < 1 || y0 > -1) { // Hit side of cylinder			
	// 		if (ray.intersection.none || t0 < ray.intersection.t_value) {
	// 			std::cout<< "yas3\n";
	// 			Point3D intersection_point = origin + t0 * dir;
	// 			ray.intersection.t_value = t0;
	// 			ray.intersection.point =  modelToWorld * intersection_point;
	// 			ray.intersection.normal = worldToModel.transpose() * Vector3D(intersection_point[0], 0, intersection_point[2]);
	// 			ray.intersection.normal.normalize();
	// 			ray.intersection.none = false;
	// 			ray.intersection.shape = "cylinder";
	// 			return true;
	// 		}
	// 	}
	// }

	// return false;
		// Unit cylinder: circle of radius 1 along the y axis between 0.5 and -0.5

	// Transform ray into object space
	Ray3D object_space_ray = Ray3D(worldToModel * ray.origin, worldToModel * ray.dir);

	double a = pow(object_space_ray.dir[0],2) + pow(object_space_ray.dir[2],2);
	double b = 2*object_space_ray.dir[0]*object_space_ray.origin[0] + 2 * object_space_ray.dir[2] * object_space_ray.origin[2];
	double c = pow(object_space_ray.origin[0],2) + pow(object_space_ray.origin[2],2) - 1;

	double k2 = b*b - 4*a*c;

	if (k2 >= 0) {
		double t1 = ((-1)*b + sqrt(k2))/(2*a);
		double t2 = ((-1)*b - sqrt(k2))/(2*a);

		double t = -1;

		// Take the lower of the two roots if both are positive (t2)
		// Take the only positive root (t1)
		if (t1 >= 0 && t2 < 0) {
			t = t1;
		} else if (t2 >= 0) {
			t = t2;
		}

		if (t != -1) {
			Point3D intersection_point = object_space_ray.origin + t*object_space_ray.dir;

			if (intersection_point[1]<=0.5 && intersection_point[1]>=-0.5) {
				if (ray.intersection.none || t < ray.intersection.t_value) {
					ray.intersection.t_value = t;
					ray.intersection.point = modelToWorld * intersection_point;
					ray.intersection.normal = worldToModel.transpose() * Vector3D(intersection_point[0], 0, intersection_point[2]);
					ray.intersection.normal.normalize();
					if (ray.dir.dot(ray.intersection.normal) > 0) {
						ray.intersection.normal = -1 * ray.intersection.normal;
					}
					ray.intersection.none = false;
					ray.intersection.shape = "cylinder";
					return true;
				}
			}
		}

		// Check intersection with the circular planes
		double t3 = (0.5 - object_space_ray.origin[1])/object_space_ray.dir[1];
		double t4 = (-0.5 - object_space_ray.origin[1])/object_space_ray.dir[1];


		if (t3 >= 0 ) {
			Point3D intersection_point = object_space_ray.origin + t3*object_space_ray.dir;
			if ((pow(intersection_point[0],2) + pow(intersection_point[2],2)) <= 1) {
				if (ray.intersection.none || t3 < ray.intersection.t_value) {
					ray.intersection.t_value = t3;
					ray.intersection.point =  modelToWorld * intersection_point;
					ray.intersection.normal = worldToModel.transpose() * Vector3D(0, 1, 0);
					ray.intersection.normal.normalize();
					if (ray.dir.dot(ray.intersection.normal) > 0) {
						ray.intersection.normal = -1 * ray.intersection.normal;
					}
					ray.intersection.none = false;
					ray.intersection.shape = "cylinder";
					return true;
				}
			}
		}

		if (t4 >= 0 ) {
			Point3D intersection_point = object_space_ray.origin + t4*object_space_ray.dir;
			if ((pow(intersection_point[0],2) + pow(intersection_point[2],2)) <= 1) {
				if (ray.intersection.none || t4 < ray.intersection.t_value) {
					ray.intersection.t_value = t4;
					ray.intersection.point =  modelToWorld * intersection_point;
					ray.intersection.normal = worldToModel.transpose() * Vector3D(0, -1, 0);
					ray.intersection.normal.normalize();
					if (ray.dir.dot(ray.intersection.normal) > 0) {
						ray.intersection.normal = -1 * ray.intersection.normal;
					}
					ray.intersection.none = false;
					ray.intersection.shape = "cylinder";
					return true;
				}
			}
		}
	}

	return false;
}