#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>


bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assert the correct values of the W coords of the origin and direction.
    // You can comment this out if you take great care to construct the rays
    // properly.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
	//!!! USEFUL NOTES: to calculate depth in localIntersect(), if the intersection happens at
	//ray.origin + ray.direction * t, then t is just the depth
	//!!! USEFUL NOTES: Here direction might be scaled, so you mustn't renormalize it in
	//localIntersect(), or you will get a depth in local coordinate system,
	//which can't be compared with intersections with other objects
    if (localIntersect(local_ray, hit)) 
	{	
		// Assert correct values of W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transform intersection coordinates into global coordinates.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}


bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
    // @@@@@@ YOUR CODE HERE
	// You could also use pure geometry relations rather than 
	//analytical tools shown in the slides
	// Here in local coordinate system, the sphere is centered on (0, 0, 0)
	//with radius 1.0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
    double px = ray.origin[0];
    double py = ray.origin[1];
    double pz = ray.origin[2];
    
    

    double vx = ray.direction[0];
    double vy = ray.direction[1];
    double vz = ray.direction[2];
    
    //printf("wait but why %f %f %f \n", vx, vy, vz);
    
    double tval;
    
    double aCoeff = pow(vx, 2) + pow(vy, 2) + pow(vz, 2);
    double bCoeff = 2 * ((px * vx)+ (py * vy) + (pz * vz));
    double cCoeff = (pow(px, 2) + pow(py, 2) + pow(pz, 2)) - (pow(this->radius, 2));
    
    double discrmnt = (pow (bCoeff, 2)) - (4 * aCoeff * cCoeff);
    bool inSphere = false;
    
    if (discrmnt < 0) {
        return false;
    } else if (discrmnt > 0) { //two solutions
        double sol1 = (- bCoeff + sqrt(discrmnt)) / (2 * aCoeff);
        double sol2 = (- bCoeff - sqrt(discrmnt)) / (2 * aCoeff);
		
		//printf("d1 %f; d2 %f\n", sol1, sol2);
        if (sol1 >= 0 && sol2 >= 0) {
			
            if (sol1 < sol2) {
                tval = sol1;
            } else {
                tval = sol2;
            }
        } else if (sol1 >= 0 && sol2 < 0) {
            inSphere = true;
            tval = sol1;
        } else if (sol2 >= 0 && sol1 < 0) {
            inSphere = true;
            tval = sol2;
        } else {
            return false;
        }
    } else if (discrmnt == 0) { // one solution
        double sol = (-bCoeff) / (2 * aCoeff);
        if (sol >= 0) {
            tval = sol;
        } else {
            return false;
        }
    }
    if (hit.depth > tval) {
        hit.depth = tval;
        
        Vector position = Vector((px + vx * tval), (py + vy * tval), (pz + vz * tval));
        hit.position = position;
        
        if (inSphere) {
            hit.normal = -position.normalized();
        } else {
            hit.normal = position.normalized();
        }
        return true;
    } else {
        return false;
    }
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	// @@@@@@ YOUR CODE HERE
	// Do not accept intersection while ray is inside the plane
	// Here in local coordinate system, the plane is at z = 0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	
	
	double px = ray.origin[0];
    double py = ray.origin[1];
    double pz = ray.origin[2];
	
    double vx = ray.direction[0];
    double vy = ray.direction[1];
    double vz = ray.direction[2];
	
	// no d, equation of plane is literally z = 0
	double tval = (-pz) / (vz);
	//printf("%f", tval);
	if (tval < 0) {
		return false;
	}
    
    if (hit.depth > tval) {
        hit.depth = tval;
	
        Vector position = Vector((px + vx * tval), (py + vy * tval), (pz + vz * tval));
        hit.position = position;
        hit.normal = Vector(0,0,1,0); //plane points up in z direction
        return true;
	} else {
        return false;
    }
}


bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
    // @@@@@@ YOUR CODE HERE (creative license)
    return false;
}


// Intersections!
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Bounding box check
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Ray parallel to bounding box plane and outside of box!
				return false;
			}
			// Ray parallel to bounding box plane and inside box: continue;
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Ensure t1 <= t2

			if (t1 > tNear) tNear = t1; // We want the furthest tNear
			if (t2 < tFar) tFar = t2; // We want the closest tFar

			if (tNear > tFar) return false; // Ray misses the bounding box.
			if (tFar < 0) return false; // Bounding box is behind the ray.
		}
	}
	// If we made it this far, the ray does intersect the bounding box.

	// The ray hits the bounding box, so check each triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}

bool Mesh::intersectTriangle(Ray const &ray,
	Triangle const &tri,
	Intersection &hit) const
{
	// Extract vertex positions from the mesh data.
	Vector const &p0 = positions[tri[0].pi];
	Vector const &p1 = positions[tri[1].pi];
	Vector const &p2 = positions[tri[2].pi];

	// @@@@@@ YOUR CODE HERE
	// Decide whether ray intersects the triangle (p0,p1,p2).
	// If so, fill in intersection information in hit and return true.
	// You may find it useful to use the routine implicitLineEquation()
	// to compute the result of the implicit line equation in 2D.
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	//!!! USEFUL NOTES: for the intersection point, its normal should satisfy hit.normal.dot(ray.direction) < 0
	
    
    double px = ray.origin[0];
    double py = ray.origin[1];
    double pz = ray.origin[2];
    double vx = ray.direction[0];
    double vy = ray.direction[1];
    double vz = ray.direction[2];
	
	//first, get normal of the triangle
	Vector e0a = p1 - p0; //edges
	Vector e0b = p2 - p0;
	Vector norm1 = e0a.cross(e0b);
    Vector norm2 = e0b.cross(e0a);
	norm1 = norm1.normalized();
    norm2 = norm2.normalized();
	Vector finalNormal;

    // if (norm1.dot(ray.direction) < 0.0) {
        // finalNormal = norm1;
    // } else if (norm2.dot(ray.direction) > 0.0) {
        // finalNormal = norm1;
    // } else {
        // return false;
    // }
    

    float dVar = norm1.dot(p0); // variable D from the plane equation in normals
    float tval = (((-norm1.dot(ray.origin)) + dVar) / norm1.dot(ray.direction));
    
    if (tval < 0) { // the ray does not intersect the plane
        return false;
    }
    
    Vector position = Vector((px + vx * tval), (py + vy * tval), (pz + vz * tval));
        
    Vector c0 = position - p0;
    Vector c1 = position - p1;
    Vector c2 = position - p2;
    
    Vector e0 = p1 - p0;
    Vector e1 = p2 - p1;
    Vector e2 = p0 - p2;
    
    if (norm1.dot(e0.cross(c0)) < 0) {
        return false;
    } 
    
    if (norm1.dot(e1.cross(c1)) < 0) {
        return false;
    }
    
    if (norm1.dot(e2.cross(c2)) < 0) {
        return false;
    }
    
    // so now we know the point is inside the triangle.
    // time to calculate the normal
    finalNormal = norm1;
        
    //if (finalNormal.dot(ray.direction) <= 0.0) {
    //    finalNormal = -finalNormal;
    //}
    
    if (hit.depth > tval) {
        //printf("%f; ", finalNormal[0]);
        hit.depth = tval;
        // point is inside triangle.
        hit.depth = tval;
        hit.position = position;
        hit.normal = finalNormal.normalized(); 
        
        return true;
	} else {
        return false;
    }
}
