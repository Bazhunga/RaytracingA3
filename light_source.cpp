/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

Colour PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	// Remember that every vector needs to be normalized. 
	// When we discuss diffusivity and specularity, we care about the angle
	// at which the light strikes the surface. We also have the values of
	// dot products be < 1, which prevents our colours from maxing out.
	// If we do not normalize, we take magnitudes into account as well, which is incorrect,
	// since we max out the colours.

	//Get the normal
	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	// Get the ray from intersection point to light
	Point3D light_pos = get_position();
	Vector3D light_ray = light_pos - ray.intersection.point;
	light_ray.normalize();

	// Get the reflected ray following formula 2*(n.s)*n - s
	double n_dot_light = normal.dot(light_ray);
	Vector3D reflected_ray = 2 * n_dot_light * normal - light_ray;
	reflected_ray.normalize();
	
	double specular_dot_product = reflected_ray.dot((-1)*ray.dir);
	double ambient_intensity = 1.0;
	double diffuse_intensity = 1.0;
	double specular_intensity = 1.2;
	double spec_shininess = ray.intersection.mat->specular_exp;

	// Need to use std::max to make sure we never have a negative colour value. 
	// Our surface is characterized by the material's diffusivity and specularity components
	Colour light_ambient = 	ambient_intensity * ray.intersection.mat->ambient;
	Colour light_diffuse = std::max(0.0, n_dot_light)*diffuse_intensity*ray.intersection.mat->diffuse;
	Colour light_specular = pow(std::max(0.0, specular_dot_product),spec_shininess)*specular_intensity*ray.intersection.mat->specular;
	// Colour material_ambient 
	// Colour material_diffuse
	// Colour material_specular	

	// ray.col = light_ambient;
	// if (n_dot_light < 0) {
	// 	return Colour(1,1,1);
	// }

	// ray.col = light_ambient + light_diffuse + light_specular;
	// ray.col.clamp();

	// ray.col = ray.col + 

	return (light_ambient + light_diffuse + light_specular);
}

