/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}
void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    traverseScene(node,ray,_modelToWorld,_worldToModel);
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, const Matrix4x4& modelToWorld, const Matrix4x4& worldToModel ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	Matrix4x4 myModelToWorld = modelToWorld*node->trans;
	Matrix4x4 myWorldToModel = node->invtrans*worldToModel;
	if (node->obj) {
		// Perform intersection.
		// Note that we pass the ray in WORLD COORDINATES at the moment
		if (node->obj->intersect(ray, myWorldToModel, myModelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray, myModelToWorld,myWorldToModel);
		childPtr = childPtr->next;
	}

}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows here if needed.
		// We know where the intersection points are
		// Now we need to draw a ray from light source to the intersection point
		// Traverse the scene
		// Now we know if that ray intersects anything
		// If it does, then we know to shade the ray darker

		Vector3D lti = curLight->light->get_position() - ray.intersection.point;
		lti.normalize();
		Point3D ip = ray.intersection.point + 0.001 * lti;
		Ray3D lightToIntersection = Ray3D(ip, lti);
		traverseScene(_root, lightToIntersection);
		if(lightToIntersection.intersection.none == false && lightToIntersection.intersection.shape.compare(ray.intersection.shape) != 0) {
			curLight->light->shade(ray, 0.8);
		}
		else {
			curLight->light->shade(ray, 1.0);
		}

		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

// The ray is in world coordinates right now
Colour Raytracer::shadeRay( Ray3D& ray ) {
	Colour currCol(0.0, 0.0, 0.0); 
	Colour refCol(0.0, 0.0, 0.0);
	Colour totalCol(0.0, 0.0, 0.0);
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	Vector3D incident = -1 * ray.dir; // We define incident as pointing away from intersection
	double selfIntersect = ray.intersection.normal.dot(incident); // negative if self intersecting
	if (!ray.intersection.none && ray.reflectNum <= _reflecNum && selfIntersect > 0) {
		computeShading(ray); 
		currCol = ray.col;  

		// Calculate the reflected ray 
		// reflectedray = 2 * (normal . incident) * normal - incident
		incident.normalize();
		// std::cout << "direction: " << ray.intersection.normal.dot(incident) << "\n";
		Vector3D reflected_vector = 2 * (ray.intersection.normal.dot(incident)) * ray.intersection.normal - incident;
		reflected_vector.normalize();
		Ray3D reflectedRay = Ray3D(ray.intersection.point + 0.001 * reflected_vector, 
													reflected_vector, 
													ray.reflectNum + 1, 
													ray.intersection.shape);
		// std::cout << "Coeff: " << ray.intersection.mat->reflection_coeff << "\n";
		// std::cout << "RN:" << ray.reflectNum  << "\n";
		// std::cout << "Power:" << pow(ray.intersection.mat->reflection_coeff, ray.reflectNum)  << "\n\n";
		refCol = shadeRay(reflectedRay);
		std::cout << "Curr Int" << ray.reflectNum << ": " << ray.intersection.shape << "\n";
		std::cout << "Next Int: " << reflectedRay.intersection.shape << "\n";
		if(reflectedRay.intersection.none || reflectedRay.intersection.shape.compare(ray.intersection.shape) == 0) { // hit nothing
			totalCol = currCol; // Only colour of the material that ray hit; no reflection
		}
		else {
			// totalCol = currCol + refCol; // Display entirely the reflected stuff
			refCol = shadeRay(reflectedRay);
			totalCol = (0.4*(currCol) + 0.6*refCol);
		}
		// std::cout << "omg: " << refCol << "\n";
	}
	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  
	totalCol.clamp();
	return totalCol; 
}	

Colour Raytracer::shadeRay2( Ray3D& ray ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray); 
		col = ray.col;  
	}
	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	_reflecNum = 4;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);


	// Sets up ray origin and direction in view space, 
	// image plane is at z = -1.
	Point3D origin(0, 0, 0);
	Point3D imagePlane;
	imagePlane[2] = -1;

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Shoot 25 instead
			// Add it to a colour object. Divide by number of rays
			Colour col_tot; 
			for (int k = -2; k <= 2; k++) {
				for (int l = -2; l<= 2; l++) {
					// This gives us the x and y coordinates of the pixel we're shooting through
					double jitterDirectionX = ((double) rand() / (RAND_MAX))/2.0;
					double jitterDirectionY = ((double) rand() / (RAND_MAX))/2.0;

					// Get rid of the 0.5 because now we use jitterDirection
					imagePlane[0] = (-double(width)/2 + j + jitterDirectionX)/factor;
					imagePlane[1] = (-double(height)/2 + i + jitterDirectionY)/factor;

					Ray3D ray;
					ray.origin = viewToWorld * origin;
					ray.dir = viewToWorld * (imagePlane - origin);//Imageplane and origin?
					ray.dir.normalize();
					col_tot = col_tot + shadeRay(ray);
				}
			}
			Colour col = (1.0/25.0) * col_tot; // Ray is in world coordinates right now 
			_rbuffer[i*width+j] = int(col[0]*255);
			_gbuffer[i*width+j] = int(col[1]*255);
			_bbuffer[i*width+j] = int(col[2]*255);
		}
	}

	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	// int width = 320; //600; //320; 
	// int height = 240;//450; //240; 
	int width = 600; //320; 
	int height = 450; //240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), 
						Colour(0.75164, 0.60648, 0.22648), 
						Colour(0.628281, 0.555802, 0.366065), 
						51.2, 0.4);

	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
						Colour(0.316228, 0.316228, 0.316228), 
						12.8, 0.0);

	Material slate( Colour(0.3, 0.3, 0.3),
						 Colour(0.0196078431, 103.0/255.0, 135.0/255.0), 
						 Colour(0.628281, 0.555802, 0.366065),
						 69, 0.5);

	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &slate );
	// SceneDagNode* planeTexture = raytracer.addObject( new UnitSquareTextured(), &slate );

	// SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &jade);
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -6));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	// raytracer.trans/late(sphere2, Vector3D(0, 1, -5));

	// >> Non textured
	raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// >> Textured
	// raytracer.translate(planeTexture, Vector3D(0, 0, -7));	
	// raytracer.rotate(planeTexture, 'z', 45); 
	// raytracer.scale(planeTexture, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	// raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	// Point3D eye2(1,0,-3);
	// Vector3D view2(-1, 0, 1);
	// Default
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);

	// Side profile
	// Point3D eye2(0, -4,2);
	// Vector3D view2(0.2, 0.4, -1);

	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	return 0;
}

