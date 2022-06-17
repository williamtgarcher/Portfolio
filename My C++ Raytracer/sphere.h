#ifndef SPHEREH
#define SPHEREH

#include "hitable.h"
#include "material.h"

void get_sphere_uv(const vec3& p, float& u, float& v) {
	float phi = atan2(p.z(), p.x());
	float theta = asin(p.y());
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / (M_PI);
}

class sphere : public hitable {
public:
	vec3 center;
	float radius;
	material* mat_ptr;

	sphere() {}
	sphere(vec3 cen, float r, material* m) : center(cen), radius(r), mat_ptr(m) {};

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	vec3 oc = r.origin - center;

	float a = r.direction.squared_length();
	float b = dot(oc, r.direction);
	float c = oc.squared_length() - radius * radius;

	// As sphere::hit is one of the most called functions, the calls to ray.origin() and ray.direction() take up a significant amount of time. 
	// Replacing these with just accessing the ray's A and B attributes directly leads to a render time reduction of ~25%
	/*vec3 oc = r.A - center;

	float a = r.B.squared_length();
	float b = dot(oc, r.B);
	float c = oc.squared_length() - radius * radius;*/

	float discriminant = b * b - a * c;
	if (discriminant > 0) {

		float temp = (-b - sqrt(discriminant)) / a; // lower root
		if (t_min < temp && temp < t_max) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			get_sphere_uv(rec.normal, rec.u, rec.v);
			rec.matptr_index = mat_ptr->id;
			return true;
		}
		//return false;///debugging only

		temp = (-b + sqrt(discriminant)) / a; // higher root
		if (t_min < temp && temp < t_max) {
			rec.t = temp;
			rec.p = r.point_at_parameter(rec.t);
			rec.normal = (rec.p - center) / radius;
			get_sphere_uv(rec.normal, rec.u, rec.v);
			rec.matptr_index = mat_ptr->id;
			return true;
		}
	}

	return false;
}


bool sphere::bounding_box(float t0, float t1, aabb& box) const {
	box = aabb(center - vec3(radius, radius, radius), center + vec3(radius, radius, radius));
	return true;
}

#endif