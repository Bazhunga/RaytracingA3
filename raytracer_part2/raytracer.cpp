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
	// >> New
	Colour totalCol = Colour(0.0,0.0,0.0);
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;

		// std::cout << "Pos light " << curLight->light->get_id() << " X: " << curLight->light->get_position()[0] << " Y: " << curLight->light->get_position()[1] <<  " Z: " << curLight->light->get_position()[2] << "\n"; 
		Vector3D lti = curLight->light->get_position() - ray.intersection.point;
		lti.normalize();
		Point3D ip = ray.intersection.point + 0.001 * lti;
		Ray3D lightToIntersection = Ray3D(ip, lti);
		traverseScene(_root, lightToIntersection);
		if(lightToIntersection.intersection.none == false && lightToIntersection.intersection.shape.compare(ray.intersection.shape) != 0) {

			curLight->light->shade(ray, 1.0);
		}
		else {
			// std::cout << "light: " << curLight->light->get_id() << "\n";
			curLight->light->shade(ray, 1.0);
		}
		totalCol = totalCol + ray.col;
		curLight = curLight->next;
	}

	if (_numlights != 0) {
		ray.col = (1.0/_numlights) * totalCol; // Get the average colour 
	} 
	else {
		ray.col = totalCol;
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
	Colour totalRef(0.0,0.0,0.0);
	Colour glossCol(0.0, 0.0, 0.0);
	Colour totalGloss(0.0,0.0,0.0);
	Colour totalCol(0.0, 0.0, 0.0);
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	Vector3D incident = -1 * ray.dir; // We define incident as pointing away from intersection
	incident.normalize();
	double selfIntersect = ray.intersection.normal.dot(incident); // negative if self intersecting
	if (!ray.intersection.none && selfIntersect > 0) {
		computeShading(ray); 
		currCol = ray.col;  

		if (ray.reflectNum <= _reflecNum && _reflectionSwitch) {
			// Calculate the reflected ray 
			// reflectedray = 2 * (normal . incident) * normal - incident
			// std::cout << "direction: " << ray.intersection.normal.dot(incident) << "\n";

			// Ensure that the incident and normal are pointing in the same direction; 
			Vector3D reflected_vector = 2 * (ray.intersection.normal.dot(incident)) * ray.intersection.normal - incident;
			reflected_vector.normalize();
			Point3D curr_int = ray.intersection.point + 0.001 * reflected_vector;
			Ray3D reflectedRay = Ray3D(curr_int, 
														reflected_vector, 
														ray.reflectNum + 1, 
														ray.intersection.shape);
			// std::cout << "Coeff: " << ray.intersection.mat->reflection_coeff << "\n";
			// std::cout << "RN:" << ray.reflectNum  << "\n";
			// std::cout << "Power:" << pow(ray.intersection.mat->reflection_coeff, ray.reflectNum)  << "\n\n";
			refCol = shadeRay(reflectedRay);
			// std::cout << "Curr Int" << ray.reflectNum << ": " << ray.intersection.shape << "\n";
			// std::cout << "Next Int: " << reflectedRay.intersection.shape << "\n";
			if(reflectedRay.intersection.none) { // hit nothing
				// Calculate how far the the reflected ray was from the light source
				// We want the reflected rays that are closer to the light source to be shaded more "phong"ly
				// The ones that fire off into the distance should be dark
				LightListNode* curLight = _lightSource;
				Point3D clp = curLight -> light -> get_position();
				Vector3D i2l = clp - curr_int;
				i2l.normalize();
				double dark_factor = reflected_vector.dot(i2l); // High if angle is small. Small if angle is large.
				totalRef = dark_factor*currCol; // Only colour of the material that ray hit; no reflection
			}
			else {
				// totalRef = currCol + refCol; // Display entirely the reflected stuff
				// refCol = shadeRay(reflectedRay);
				// totalRef = (0.4*(currCol) + 0.6*refCol);
				totalRef = (0.4*currCol + 0.8* pow(ray.intersection.mat->reflection_coeff, ray.reflectNum) * refCol);
				// totalRef = 0.2*currCol + 0.9*refCol;
			}
		}
		if (ray.reflectNum <= _reflecNum && _glossSwitch) {
			Vector3D reflected_vector = 2 * (ray.intersection.normal.dot(incident)) * ray.intersection.normal - incident;
			reflected_vector.normalize();
			Point3D curr_int = ray.intersection.point + 0.001 * reflected_vector;

			Point3D ip = ray.intersection.point + 0.001 * reflected_vector;
			// Generate a random vector
			double randX = ((double) rand() / (RAND_MAX));
			double randY = ((double) rand() / (RAND_MAX));
			double randZ = ((double) rand() / (RAND_MAX));
			Vector3D randoVec = Vector3D(randX, randY, randZ);

			Vector3D pv_u = reflected_vector.cross(randoVec);
			Vector3D pv_v = reflected_vector.cross(pv_u);
			pv_u.normalize();
			pv_v.normalize();

			// Now we have the equations of the plane that's perpendicular to the reflected_vector
			// We want to get the shading from all of these vectors and then average them
			int a = 9;
			for (int i = 0; i < a; i++ ) {
				double jitterU = ((double) rand() / (RAND_MAX)); // from 0 to a
				double jitterV = ((double) rand() / (RAND_MAX)); // from 0 to a
				if (((double) rand() / (RAND_MAX)) > 0.5){ jitterU = -1 * jitterU;}
				if (((double) rand() / (RAND_MAX)) > 0.5){ jitterV = -1 * jitterV;}
				// std::cout<<"ju " << jitterU << "\n";
				// std::cout<<"jv " << jitterV << "\n";
				Vector3D jitteredVector = 2.8* reflected_vector + jitterU*pv_u + jitterV*pv_v;
				jitteredVector.normalize();

				Ray3D glossRay = Ray3D(ip, 
											jitteredVector, 
											ray.reflectNum + 1, 
											ray.intersection.shape);
				glossCol = shadeRay(glossRay);
				if(glossRay.intersection.none) {
					LightListNode* curLight = _lightSource;
					Point3D clp = curLight -> light -> get_position();
					Vector3D i2l = clp - curr_int;
					i2l.normalize();
					double dark_factor = 0.5 + 0.5 * reflected_vector.dot(i2l)*ray.intersection.mat->gloss_coeff;
					totalGloss = totalGloss + dark_factor * currCol; // Only colour of the material that ray hit; no reflection
				}
				else {
					// std::cout << "r: " << glossCol[0] << "g: " << glossCol[1] << "b: " << glossCol[2]  < "\n";
					totalGloss = totalGloss + (0.8*currCol + 0.4* pow(ray.intersection.mat->gloss_coeff, ray.reflectNum) * glossCol);
					// totalGloss = totalGloss + currCol;
				}
			}
			totalGloss = (1/double(a)) * totalGloss;
		}
		if(!_reflectionSwitch && !_glossSwitch) {
			totalCol = currCol;
		}
		else {
			double tone_down = 1;
			if (_reflectionSwitch && _glossSwitch) {
				tone_down = 0.6;
			}
			totalCol = tone_down * (totalRef + totalGloss);
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
		Vector3D up, double fov, char* fileName, int numLights) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	_reflecNum = 2;
	_numlights = numLights;
	_reflectionSwitch = true;
	_glossSwitch = true;
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
		std::cout << "Drawing row: " << i << "\n";
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
	bool _shadowSwitch = false;
	// int width = 320; //600; //320; 
	// int height = 240;//450; //240; 
	int width = 320; //320; 
	int height = 150; //240; 

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
						70, 0.5, 0.8);

	Material jade( Colour(0.3, 0.3, 0.3), Colour(0.54, 0.89, 0.63), 
						Colour(0.316228, 0.316228, 0.316228), 
						12.8, 0.6, 0.8);

	Material skyblue( Colour(0.3, 0.3, 0.3),
						 Colour(0.0196078431, 103.0/255.0, 135.0/255.0), 
						 Colour(0.316228, 0.316228, 0.316228), 
						 // Colour(0.628281, 0.555802, 0.366065),
						 69, 0.8, 0.2);

	Material slate( Colour(0.3, 0.3, 0.3),
						 Colour(62.0/255, 89.0/255.0, 89.0/255.0), 
						 Colour(0.7, 0.7, 0.7), 
						 69, 0.8, 0.2);

	Material frosty( Colour (0.3, 0.3, 0.3), 
						 Colour(123.0/255.0, 136.0/255.0, 140.0/255.0),
						 Colour(0.628281, 0.555802, 0.366065),
						 90, 0.6, 0.8);

	// Defines a point light source.
	double main_light_X = 0;
	double main_light_Y = 0;
	double main_light_Z = 3;
	double toRadian = 2*M_PI/360.0;
	raytracer.addLightSource( new PointLight(Point3D(main_light_X, main_light_Y, main_light_Z), Colour(0.9, 0.9, 0.9), 1));
	// std::cout << "Coords: " << main_light_X << " " << main_light_Y << " " << main_light_Z << "\n";
	// Define an area light source as a collection of point lights
	int light_ring_num = 0;
	int outer_ring_num = 0;
	if (_shadowSwitch == true) {
		double radius = 1.0; // From center point light to surrounding point lights
		int light_ring_num = 50; // Number of lights in a ring around the center light
		double angleSpacing = 360 / double(light_ring_num); // angle spacing between lights (5 lights would have 72 degrees between each light)
		double currAngle = 0;
		for (int i = 0; i < light_ring_num; i++) {
			double currAngle = i * angleSpacing;
			double sublight_X = radius * cos(currAngle * toRadian);
			double sublight_Y = radius * sin(currAngle * toRadian); 
			// std::cout << "ILRCoords: " << sublight_X << ", " << sublight_Y << ", " << main_light_Z << "\n";
			raytracer.addLightSource( new PointLight(Point3D(sublight_X, sublight_Y, main_light_Z), 
					Colour(0.7,0.7,0.7), i + 2) );
		}

		// Define moar points
		radius = 1.5;
		int outer_ring_num = 100;
		angleSpacing = 360 / double(outer_ring_num); // angle spacing between lights (5 lights would have 72 degrees between each light)
		for (int i = 0; i < outer_ring_num; i++) {
			double currAngle = i * angleSpacing;
			double sublight_X = radius * cos(currAngle * toRadian);
			double sublight_Y = radius * sin(currAngle * toRadian); 
			// std::cout << "ORCoords: " << sublight_X << ", " << sublight_Y << ", " << main_light_Z << "\n";
			raytracer.addLightSource( new PointLight(Point3D(sublight_X, sublight_Y, main_light_Z), 
					Colour(0.7,0.7,0.7), i + light_ring_num + 2) );
		}
	}


	// Add a unit square into the scene with material mat.
	// SceneDagNode* plane2 = raytracer.addObject( new UnitSquare2(), &jade );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	// SceneDagNode* cylinder = raytracer.addObject( new UnitCylinder(), &gold );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere2(), &frosty);
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &skyblue );
	// SceneDagNode* planeTexture = raytracer.addObject( new UnitSquareTextured(), &skyblue );

	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor3[3] = { 0.5, 0.5, 0.5 };
	double factor2[3] = { 6.0, 6.0, 6.0 };

	// >> Draw Ellipsoid
	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	// >> Draw Spheres
	raytracer.translate(sphere2, Vector3D(0, 3, -5));	

	// >> Draw Cylinder
	// raytracer.translate(cylinder, Vector3D(0, 0, -4.5));	
	// // raytracer.rotate(cylinder, 'x', -45); 
	// // raytracer.rotate(cylinder, 'z', 45); 
	// raytracer.scale(cylinder, Point3D(0, 0, 0), factor1);

	// >> Draw plane
	// raytracer.translate(plane2, Vector3D(0, 0, -3));	
	// raytracer.scale(plane2, Point3D(0, 0, 0), factor3);

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

	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp", light_ring_num + outer_ring_num);
	
	return 0;
}

