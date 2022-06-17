#pragma once

#include "hitable.h"
#include "material.h"
#include "texture.h"


class isotropic : public material {
public:
	texture* albedo;

	isotropic(texture* a) : albedo(a) {}
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
		scattered = ray(rec.p, random_in_unit_sphere());
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}
};

//For accurate Rayleigh scattering, blue light must have a shorter mean free path

//class rayleigh : public material {
//public:
//	vec3 rgb_factors = vec3(0.29679, 0.58467, 1.0000);
//
//	rayleigh() {}
//	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
//		scattered = ray(rec.p, random_in_unit_sphere());
//		float intensity = 0.5 * (1 + pow(dot(scattered.direction, r_in.direction), 2));
//		attenuation = intensity * rgb_factors;
//		return true;
//	}
//};


class constant_medium : public hitable {
public:
	hitable* boundary;
	float density;
	material* phase_function;

	constant_medium(hitable* b, float d, texture* a) : boundary(b), density(d) { phase_function = new isotropic(a); }
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		return boundary->bounding_box(t0, t1, box);
	}
};

bool constant_medium::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	hit_record rec1, rec2;
	if (boundary->hit(r, -FLT_MAX, FLT_MAX, rec1)) {
		if (boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2)) {
			rec1.t = (rec1.t < t_min) ? t_min : rec1.t;
			rec2.t = (rec2.t > t_max) ? t_max : rec2.t;
			if (rec1.t >= rec2.t) {
				return false;
			}
			rec1.t = (rec1.t < 0) ? 0 : rec1.t;

			float r_direction_length = r.direction.length();
			float distance_inside_boundary = (rec2.t - rec1.t) * r_direction_length;
			float hit_distance = -(1 / density) * log(random());
			if (hit_distance < distance_inside_boundary) {
				rec.t = rec1.t + hit_distance / r_direction_length;
				rec.p = r.point_at_parameter(rec.t);
				rec.normal = vec3(1, 0, 0); // (arbitrary!)
				rec.matptr_index = phase_function->id;
				return true;
			}
		}
	}
	return false;
}