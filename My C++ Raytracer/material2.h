#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

#include "material.h"
#include "cone_direction.h"

float projectedAreaP(vec3 w_i, vec3 w_g, vec3 w_p) {
	return dot(w_i, w_p) / dot(w_p, w_g);
}
float projectedAreaT(vec3 w_i, vec3 w_g, vec3 w_p, vec3 w_t) {
	float wp_dot_wg = dot(w_p, w_g);
	return dot(w_i, w_t) * sqrt(1 - min(wp_dot_wg * wp_dot_wg, 1.0f)) / wp_dot_wg;
}


float G_1(vec3 w_i, vec3 w_m, vec3 w_g, float total_area) {
	return float(dot(w_i, w_m) > 0) * min(1.0f, dot(w_i, w_g) / total_area);
}


vec3 material2Reflection(vec3 w_r, vec3 local_normal, vec3 mer_value, vec3& attenuation, vec3 original_attenuation) {
	vec3 reflected = reflect(w_r, local_normal);
	if (mer_value[2] > 0.0) {
		if (random() < mer_value[2]) {
			attenuation *= original_attenuation;
			return local_normal + random_unit_vector();
		}
		vec3 direction;
		direction = reflected + randCone(reflected, mer_value[2] * M_PI);
		attenuation *= ray(original_attenuation, vec3(1, 1, 1) - original_attenuation).point_at_parameter(1 - mer_value[0]);
		return direction;
	}
	attenuation *= ray(original_attenuation, vec3(1, 1, 1) - original_attenuation).point_at_parameter(1 - mer_value[0]);
	return reflected;
}



vec3 microfacetAlgorithm(vec3 w_r, vec3 w_p, vec3 mer_value, const hit_record& rec, vec3& attenuation) {
	vec3 w_t = dot(w_p, rec.normal) * rec.normal - w_p;
	vec3 w_g = rec.normal;

	vec3 original_attenuation = attenuation;
	attenuation = vec3(1, 1, 1);

	float a_p = projectedAreaP(-w_r, w_g, w_p);
	float a_t = projectedAreaT(-w_r, w_g, w_p, w_t);
	float total_area = a_p + a_t;
	float prob_p = a_p / total_area;
	float prob_t = a_t / total_area;

	float ranfloat = random();
	int m = ranfloat <= prob_p; // 1:w_p, 0:w_t
	vec3 w_m;
	for (int i = 0; i < 10000; i++) {
		w_m = m ? w_p : w_t;

		w_r = material2Reflection(w_r, w_m, mer_value, attenuation, original_attenuation);
		float mask = G_1(w_r, w_m, w_g, total_area);

		if (random() <= mask) {
			return w_r;
		}
		else {
			m = 1 - m;
			a_p = projectedAreaP(-w_r, w_g, w_p);
			a_t = projectedAreaT(-w_r, w_g, w_p, w_t);
			total_area = a_p + a_t;
		}
	}
	return w_r;
}


class opaque : public material {
public:
	texture* albedo;
	texture* normal;
	texture* mer;

	opaque(texture* albedo_map, texture* normal_map, texture* mer_map) : albedo(albedo_map), normal(normal_map), mer(mer_map) {}
	virtual vec3 emitted(float u, float v, const vec3& p) const { return (pow(10, 2 * mer->value(u, v, p).e[1]) - 1) * albedo->value(u, v, p); }
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
		if (-dot(r_in.direction, rec.normal) < 0.0) {
			return false;
		}
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		vec3 normal_offset = unit_vector(2 * normal->value(rec.u, rec.v, rec.p) - vec3(1, 1, 1));
		vec3 local_normal = relativeVector(rec.normal, normal_offset);
		vec3 mer_value = mer->value(rec.u, rec.v, rec.p);
		
		vec3 direction = microfacetAlgorithm(r_in.direction, local_normal, mer_value, rec, attenuation);
		scattered = ray(rec.p, direction);
		if (attenuation.is_zero()) { return false; }
		return true;
	}
};



class translucent : public material {
public:
	texture* albedo;
	texture* normal;
	texture* roughness;
	float ref_idx;

	translucent(texture* albedo_map, texture* normal_map, texture* rough_map, float ri) : albedo(albedo_map), normal(normal_map), roughness(rough_map), ref_idx(ri) {}
	virtual vec3 emitted(float u, float v, const vec3& p) const { return (pow(10, 2 * roughness->value(u, v, p).e[1]) - 1) * albedo->value(u, v, p); }
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
		attenuation = albedo->value(rec.u, rec.v, rec.p);
		vec3 normal_offset = unit_vector(2 * normal->value(rec.u, rec.v, rec.p) - vec3(1, 1, 1));
		vec3 local_normal = relativeVector(rec.normal, normal_offset);
		vec3 roughness_value = roughness->value(rec.u, rec.v, rec.p);
		roughness_value.e[0] = 1.0;


		vec3 outward_normal;
		float ni_over_nt;
		vec3 refracted;
		float reflect_prob;
		float cosine;
		if (dot(r_in.direction, local_normal) > 0) {
			outward_normal = -local_normal;
			ni_over_nt = ref_idx;
			cosine = ref_idx * dot(r_in.direction, local_normal) / r_in.direction.length();
		}
		else {
			outward_normal = local_normal;
			ni_over_nt = 1.0 / ref_idx;
			cosine = -dot(r_in.direction, local_normal) / r_in.direction.length();
		}
		if (refract(r_in.direction, outward_normal, ni_over_nt, refracted)) {
			reflect_prob = shlick(cosine, ref_idx);
		}
		else {
			reflect_prob = 1.0;
		}

		if (random() < reflect_prob) {
			vec3 direction = microfacetAlgorithm(r_in.direction, local_normal, roughness_value, rec, attenuation);
			scattered = ray(rec.p, direction);
		}
		else {
			scattered = ray(rec.p, refracted + randCone(refracted, roughness_value[2] * M_PI));
		}

		return true;
	}
};
