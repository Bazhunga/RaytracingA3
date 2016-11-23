/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	// Remember that every vector needs to be normalized. 
	// When we discuss diffusivity and specularity, we care about the angle
	// at which the light strikes the surface. If we do not normalize, 
	// we take magnitudes into account as well, which is incorrect
	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	Point3D light_pos = get_position();
	Vector3D light_ray = light_pos - ray.intersection.point;
	light_ray.normalize();

	Vector3D reflected_ray = 2*n_dot_light*normal - light_ray;
	reflected_ray.normalize();
	
	double n_dot_light = normal.dot(light_ray);
	double specular_dot_product = reflected_ray.dot((-1)*ray.dir);
	double ambient_intensity = 1;
	double diffuse_intensity = 1;
	double specular_intensity = 1;
	double spec_shininess = ray.intersection.mat->specular_exp;

	Colour light_ambient = 	ambient_intensity* ray.intersection.mat->ambient * _col_ambient;
	Colour light_diffuse = std::max(0.0, n_dot_light)*diffuse_intensity*ray.intersection.mat->diffuse * _col_diffuse;
	Colour light_specular = pow(std::max(0.0, specular_dot_product),spec_shininess)*specular_intensity*ray.intersection.mat->specular * _col_specular;
	// Colour material_ambient 
	// Colour material_diffuse
	// Colour material_specular	

	ray.col = light_ambient;
	if (n_dot_light < 0) {
		return;
	}

	ray.col = ray.col + light_diffuse + light_specular;
	ray.col.clamp();

	// ray.col = ray.col + 

	return;
}

