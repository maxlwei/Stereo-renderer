#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#include "raytracer.hpp"
#include "image.hpp"


void Raytracer::render(const char *filename, const char *depth_filename,
                       Scene const &scene)
{
    // Allocate the two images that will ultimately be saved.
    Image colorImage(scene.resolution[0], scene.resolution[1]);
    Image depthImage(scene.resolution[0], scene.resolution[1]);
    
    // Create the zBuffer.
    double *zBuffer = new double[scene.resolution[0] * scene.resolution[1]];
    for(int i = 0; i < scene.resolution[0] * scene.resolution[1]; i++) {
        zBuffer[i] = DBL_MAX;
    }

	// @@@@@@ YOUR CODE HERE
	// calculate camera parameters for rays, refer to the slides for details
	//!!! USEFUL NOTES: tan() takes rad rather than degree, use deg2rad() to transform
	//!!! USEFUL NOTES: view plane can be anywhere, but it will be implemented differently,
	//you can find references from the course slides 22_GlobalIllum.pdf
    
    
	Vector cameraPos = scene.camera.position;
	Vector cameraCenter = scene.camera.center;
	
	Vector cameraPosR = scene.camera.position;
	Vector cameraCenterR = scene.camera.center;
	
	// viewing direction vector get by taking center and subtracting camera position
	Vector wVecOriginal = scene.camera.center; - cameraPos;
    wVecOriginal.normalize();
	// up vector is defined (u)
    Vector uVec = scene.camera.up;
    uVec.normalize();
	// right vector is gotten by taking the cross product of w and v
	Vector rVecOriginal = wVecOriginal.cross(uVec);
	rVecOriginal.normalize();
	
	double stereoDisplacement = scene.camera.stereoDist / 2.0;
	int widthResolution = scene.resolution[0];
	if (scene.camera.stereoDist > 0.0) {		
		printf("Start left picture.\n");
		cameraPos = scene.camera.position + (rVecOriginal * stereoDisplacement);
		cameraPosR = scene.camera.position - (rVecOriginal * stereoDisplacement);
		
		widthResolution = floor(scene.resolution[0] / 2);
	} else if (scene.camera.stereoDist < 0.0) {		
		printf("Start left picture.\n");
		stereoDisplacement = - scene.camera.stereoDist / 2.0;
		cameraPos = scene.camera.position - (rVecOriginal * stereoDisplacement);
		cameraPosR = scene.camera.position + (rVecOriginal * stereoDisplacement);
		
		widthResolution = floor(scene.resolution[0] / 2);
	}

	
    Vector wVec = cameraCenter - cameraPos;
    wVec.normalize();
    Vector rVec = wVec.cross(uVec);
    rVec.normalize();
    
    // get top from tan(fovy)
    double tangent = tan(deg2rad(scene.camera.fovy/2));
	//double atangent = atan(deg2rad(scene.camera.fovy)/2);
    // get length of top from centre of image plane
    double top = scene.camera.zNear * tangent;
    double right = top * scene.camera.aspect;
	if (scene.camera.stereoDist != 0.0) {		
		right = right / 2;
	}
    double left = -right;
    double bottom = -top;
	
    // calculate vector from camera to left top of image plane
    Vector centerVec = cameraPos + (scene.camera.zNear * wVec);
    Vector oVec = centerVec + (left * rVec) + (bottom * uVec);
    double deltaU = (right - left) / scene.resolution[0];
	if (scene.camera.stereoDist != 0.0) {		
		deltaU = deltaU * 2;
	}
    double deltaV = (top - bottom) / scene.resolution[1];    
	    
    // Iterate over all the pixels in the image.
    for(int y = 0; y < scene.resolution[1]; y++) {
        for(int x = 0; x < widthResolution; x++) {

            // Generate the appropriate ray for this pixel
			Ray ray;
			if (scene.objects.empty())
			{
				//no objects in the scene, then we render the default scene:
				//in the default scene, we assume the view plane is at z = 640 with width and height both 640
				ray = Ray(cameraPos, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - cameraPos).normalized());
			}
			else
			{
				// set primary ray using the camera parameters
				//!!! USEFUL NOTES: all world coordinate rays need to have a normalized direction
				
				Vector changeU = (x + 0.5) * deltaU * rVec;
                Vector changeY = (y + 0.5) * deltaV * uVec;
                Vector pixelPos = oVec + changeU + changeY;
                
                Vector rayOfHope = pixelPos - cameraPos;
                rayOfHope.normalize();
                
                ray = Ray(cameraPos, rayOfHope);
                //!!! rays do not have w coordinate constructed properly.
			}

            // Initialize recursive ray depth.
            int rayDepth = 0;
           
            // Our recursive raytrace will compute the color and the z-depth
            Vector color;

            // This should be the maximum depth, corresponding to the far plane.
            // NOTE: This assumes the ray direction is unit-length and the
            // ray origin is at the camera position.
            double depth = scene.camera.zFar;

            // Calculate the pixel value by shooting the ray into the scene
            trace(ray, rayDepth, scene, color, depth);

            // Depth test
            if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
                depth < zBuffer[x + y*scene.resolution[0]]) {
                zBuffer[x + y*scene.resolution[0]] = depth;

                // Set the image color (and depth)
                colorImage.setPixel(x, y, color);
                depthImage.setPixel(x, y, (depth-scene.camera.zNear) / 
                                        (scene.camera.zFar-scene.camera.zNear));
            }
        }

		//output step information
		if (y % 100 == 0)
		{
			printf("Row %d pixels finished.\n", y);
		}
    }
	
	if (scene.camera.stereoDist != 0.0) {		
		printf("Start right picture.\n");
		Vector wVecR = cameraCenterR - cameraPosR;
		wVecR.normalize();
		// up vector is defined (u)
		// right vector is gotten by taking the cross product of w and v
		Vector rVecR = wVecR.cross(uVec);
		rVecR.normalize();
		
		// calculate vector from camera to left top of image plane
		Vector centerVecR = cameraPosR + (scene.camera.zNear * wVecR);
		Vector oVecR = centerVecR + (left * rVecR) + (bottom * uVec);
		
		// Iterate over all the pixels in the image.
		for(int y = 0; y < scene.resolution[1]; y++) {
			for(int x = 0; x < (scene.resolution[0] / 2); x++) {

				// Generate the appropriate ray for this pixel
				Ray ray;
				if (scene.objects.empty())
				{
					//no objects in the scene, then we render the default scene:
					//in the default scene, we assume the view plane is at z = 640 with width and height both 640
					ray = Ray(cameraPosR, (Vector(-320, -320, 640) + Vector(x + 0.5, y + 0.5, 0) - cameraPosR).normalized());
				}
				else
				{
					// set primary ray using the camera parameters
					//!!! USEFUL NOTES: all world coordinate rays need to have a normalized direction
					
					Vector changeU = (x + 0.5) * deltaU * rVecR;
					Vector changeY = (y + 0.5) * deltaV * uVec;
					Vector pixelPos = oVecR + changeU + changeY;
					
					Vector rayOfHope = pixelPos - cameraPosR;
					rayOfHope.normalize();
					ray = Ray(cameraPosR, rayOfHope);
					//!!! rays do not have w coordinate constructed properly.
				}

				// Initialize recursive ray depth.
				int rayDepth = 0;
			   
				// Our recursive raytrace will compute the color and the z-depth
				Vector color;

				// This should be the maximum depth, corresponding to the far plane.
				// NOTE: This assumes the ray direction is unit-length and the
				// ray origin is at the camera position.
				double depth = scene.camera.zFar;

				// Calculate the pixel value by shooting the ray into the scene
				trace(ray, rayDepth, scene, color, depth);

				// Depth test
				int testDepth = x + floor(scene.resolution[0] / 2) + y*scene.resolution[0];
				if(depth >= scene.camera.zNear && depth <= scene.camera.zFar && 
					depth < zBuffer[testDepth]) {
					zBuffer[testDepth] = depth;

					// Set the image color (and depth)
					colorImage.setPixel(x+floor(scene.resolution[0] / 2), y, color);
					depthImage.setPixel(x+floor(scene.resolution[0] / 2), y, (depth-scene.camera.zNear) / 
											(scene.camera.zFar-scene.camera.zNear));
				}
			}

			//output step information
			if (y % 100 == 0)
			{
				printf("Row %d pixels finished.\n", y);
			}
		}
	}

	//save image
    colorImage.writeBMP(filename);
    depthImage.writeBMP(depth_filename);

	printf("Ray tracing finished with images saved.\n");

    delete[] zBuffer;
}


bool Raytracer::trace(Ray const &ray, 
                 int &rayDepth,
                 Scene const &scene,
                 Vector &outColor, double &depth)
{
    // Increment the ray depth.
    rayDepth++;

    // - iterate over all objects calling Object::intersect.
    // - don't accept intersections not closer than given depth.
    // - call Raytracer::shade with the closest intersection.
    // - return true iff the ray hits an object.
	if (scene.objects.empty())
	{
		// no objects in the scene, then we render the default scene:
		// For default, we assume there's a cube centered on (0, 0, 1280 + 160) with side length 320 facing right towards the camera
		// test intersection:
		double x = 1280 / ray.direction[2] * ray.direction[0] + ray.origin[0];
		double y = 1280 / ray.direction[2] * ray.direction[1] + ray.origin[1];
		if ((x <= 160) && (x >= -160) && (y <= 160) && (y >= -160))
		{
			//if intersected:
			Material m; m.emission = Vector(16.0, 0, 0); m.reflect = 0; //just for default material, you should use the intersected object's material
			Intersection intersection;	//just for default, you should pass the intersection found by calling Object::intersect()
			outColor = shade(ray, rayDepth, intersection, m, scene);
			depth = 1280;	//the depth should be set inside each Object::intersect()
		}
	}
	else
	{
		// Note that for Object::intersect(), the parameter hit is the current hit
		// your intersect() should be implemented to exclude intersection far away than hit.depth
		Intersection hit;
		hit.depth = depth;
		Material m;
		
        for(auto it = scene.objects.begin(); it != scene.objects.end(); ++ it) {
            Object * curObj = *it;
            
            if (curObj->intersect(ray, hit)) {
				m = curObj->material;
                outColor = shade(ray, rayDepth, hit, m, scene);
            }
        }
        
        depth = hit.depth;
	}

    // Decrement the ray depth.
    rayDepth--;

    return false; 
}


Vector Raytracer::shade(Ray const &ray,
                 int &rayDepth,
                 Intersection const &intersection,
                 Material const &material,
                 Scene const &scene)
{
    // - iterate over all lights, calculating ambient/diffuse/specular contribution
    // - use shadow rays to determine shadows
    // - integrate the contributions of each light
    // - include emission of the surface material
    // - call Raytracer::trace for reflection/refraction colors
    // Don't reflect/refract if maximum ray recursion depth has been reached!
	//!!! USEFUL NOTES: attenuate factor = 1.0 / (a0 + a1 * d + a2 * d * d)..., ambient light doesn't attenuate, nor does it affected by shadow
	//!!! USEFUL NOTES: don't accept shadow intersection far away than the light position
	//!!! USEFUL NOTES: for each kind of ray, i.e. shadow ray, reflected ray, and primary ray, the accepted furthest depth are different
	Vector diffuse(0);
	Vector ambient(0);
	Vector specular(0);
	Vector refractedLight(0);
	    
    for (auto lightIter = scene.lights.begin(); lightIter != scene.lights.end(); lightIter++)
	{
		// calculate local illumination here, remember to add all the lights together
		// also test shadow here, if a point is in shadow, multiply its diffuse and specular color by (1 - material.shadow)
		
		Vector newAmbient = (lightIter->ambient * material.ambient);
		
		Vector normalI = intersection.normal.normalized();
		Vector incomingLight = (lightIter->position - intersection.position).normalized();
				
		float diffuseDot = max(normalI.dot(incomingLight), 0.0);
		Vector newDiffuse = (diffuseDot * material.diffuse * lightIter->diffuse);
		
		
		Vector reflectedLight = (-2 * (incomingLight.dot(normalI))* normalI + incomingLight).normalized();
		Vector viewVector = ray.direction.normalized();
		
		float specularDot = max(viewVector.dot(reflectedLight), 0.0);
		Vector newSpecular = (pow(specularDot, material.shininess) * material.specular * lightIter->specular);
		
        // create shadow ray
        
		Intersection shdHit;
        Ray shadowRay = Ray((intersection.position + (incomingLight * 0.000001)), incomingLight);
        double lightDist = (lightIter->position - intersection.position).length();
        shdHit.depth = lightDist;
		bool hitShadow = false;
        
        for(auto it = scene.objects.begin(); it != scene.objects.end(); ++ it) {
            Object * curObj = *it;
            
            if (curObj->intersect(shadowRay, shdHit)) {
				hitShadow = true;
            }
        }
        
        
        double attFactor = (lightIter->attenuation[0] + lightIter->attenuation[1] * lightDist + lightIter->attenuation[2] * lightDist * lightDist);
        
        if (hitShadow) {
            newSpecular = newSpecular * (1 - material.shadow);
            newDiffuse = newDiffuse * (1 - material.shadow);
        }
        
        newSpecular = newSpecular / attFactor;
        newDiffuse = newDiffuse  / attFactor;
        
        specular += newSpecular;
        ambient += newAmbient;
        diffuse += newDiffuse;
        
		
	}

	Vector reflectedLight(0);
	if ((!(ABS_FLOAT(material.reflect) < 1e-6)) && (rayDepth < MAX_RAY_RECURSION))
	{
		if (!ray.refracting) { //inside a refracting object
            // calculate reflected color using trace() recursively
            Vector normalI = intersection.normal.normalized();
            Vector reflectedDir = (-2 * (ray.direction.dot(normalI))* normalI + ray.direction).normalized();
            Ray reflectRay = Ray((intersection.position + (reflectedDir * 0.000001)), reflectedDir);
            double superdepth = DBL_MAX;
            trace(reflectRay, rayDepth, scene, reflectedLight, superdepth);
        }
	}
    
    if ((!(ABS_FLOAT(material.refract) < 1e-6)) && (rayDepth < MAX_RAY_RECURSION)) {
        
        Vector normalI = intersection.normal.normalized();
        
        double refractn = 1.0;
        if (ray.refracting) { //inside a refracting object
            refractn = material.rfrIndex / 1.0;
        } else { 
            refractn = 1.0 / material.rfrIndex;
        }
        double cosI = normalI.dot(ray.direction);
        double sinT2 = refractn * refractn * (1.0 - cosI * cosI);
        if (sinT2 <= 1.0)
        {
            Vector refractDir = refractn * ray.direction - (refractn + sqrt(1.0 - sinT2)) * normalI;
            Ray refractRay = Ray((intersection.position + (refractDir * 0.001)), refractDir, !(ray.refracting));
            double superdepth = DBL_MAX;
            trace(refractRay, rayDepth, scene, refractedLight, superdepth);
        }
	}

	return material.emission + ambient + diffuse + specular + material.reflect * reflectedLight + material.refract * refractedLight;
}