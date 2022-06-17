#pragma once

#include "hitable.h"
#include "material.h"
#include "texture.h"

//class straight_on : public material {
//public:
//	texture* albedo;
//
//	straight_on(texture* a) : albedo(a) {}
//	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
//		scattered = ray(rec.p, r_in.direction);
//		attenuation = rec.normal * albedo->value(rec.u, rec.v, rec.p);
//		return true;
//	}
//};
//
//
//class medium : public hitable {
//public:
//	hitable* boundary;
//	float density;
//	material* phase_function;
//
//	medium(hitable* b, float d, texture* a) : boundary(b), density(d) { phase_function = new straight_on(a); }
//	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
//	virtual bool bounding_box(float t0, float t1, aabb& box) const {
//		return boundary->bounding_box(t0, t1, box);
//	}
//};
//
//bool medium::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
//	hit_record rec1, rec2;
//	if (boundary->hit(r, -FLT_MAX, FLT_MAX, rec1)) {
//		if (boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2)) {
//			rec1.t = (rec1.t < t_min) ? t_min : rec1.t;
//			rec2.t = (rec2.t > t_max) ? t_max : rec2.t;
//			if (rec1.t >= rec2.t) {
//				return false;
//			}
//			rec1.t = (rec1.t < 0) ? 0 : rec1.t;
//
//
//			float r_direction_length = r.direction.length();
//			float distance_inside_boundary = (rec2.t - rec1.t) * r_direction_length;
//			rec.t = rec2.t;
//			rec.p = r.point_at_parameter(rec.t);
//			//rec.normal = vec3(1, 0, 0); // (arbitrary!)
//			float attenuation = exp(-density * distance_inside_boundary);
//			rec.normal = attenuation * vec3(1, 1, 1); // dont need normal so using it to carry attenuation from medium
//			rec.mat_ptr = phase_function;
//			return true;
//		}
//	}
//	return false;
//}

class medium : public material {
public:
	vec3 absorption;
	vec3 emission;

	medium(vec3 albedo, float density, vec3 emittance, float brightness) : absorption(albedo * density), emission(emittance * brightness) { is_medium = true; }
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
		scattered = ray(rec.p, r_in.direction);
		attenuation = vec3(1, 1, 1);

		scattered.media_inside = r_in.media_inside;
		if (dot(r_in.direction, rec.normal) < 0) {
			scattered.media_inside[id]++;
		}
		else {
			scattered.media_inside[id]--;
		}

		return true;
	}
	
	static vec3 media_absorption(const ray& r_in, const hit_record& rec) {
		vec3 total_absorption(0, 0, 0);
		for (int i = 0; i < m; i++) {
			if (r_in.media_inside[i] > 0) {
				medium* medium_pointer = (medium*)pointer_array[i];
				total_absorption += medium_pointer->absorption;
			}
		}
		return exp(-rec.t * r_in.direction.length() * total_absorption);
	}

	static vec3 media_emission(const ray& r_in, const hit_record& rec) {
		vec3 total_emission(0, 0, 0);
		for (int i = 0; i < m; i++) {
			if (r_in.media_inside[i] > 0) {
				medium* medium_pointer = (medium*)pointer_array[i];
				total_emission += (medium_pointer->emission / medium_pointer->absorption) * (vec3(1, 1, 1) - exp(-rec.t * r_in.direction.length() * medium_pointer->absorption));
			}
		}
		return total_emission;
	}
};