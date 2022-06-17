#pragma once

#define _USE_MATH_DEFINES

#include "vec3.h"
#include "random.h"
#include <math.h>

vec3 randCone(float half_angle) {
	float z_min = cos(half_angle);
	float z = z_min + (1 - z_min) * random();
	float theta = 2 * M_PI * random();
	float r = sqrt(1 - z * z);
	float x = r * cos(theta);
	float y = r * sin(theta);
	return vec3(x, y, z);
}

vec3 randCone(vec3 direction, float half_angle) {
	vec3 a(0, 0, 1);
	vec3 b = direction / direction.length();
	//vec3 u = cross(a, b);
	vec3 u(-b[1], b[0], 0);
	float sin_phi = min(u.length(), 1.0f);
	u /= sin_phi;
	float phi = asin(sin_phi);
	float cos_phi = cos(phi);
	float not_cos_phi = 1 - cos_phi;

	vec3 randvec = randCone(half_angle);
	float x = (cos_phi + not_cos_phi * u[0] * u[0]) * randvec[0] + (not_cos_phi * u[0] * u[1]) * randvec[1] + (sin_phi * u[1]) * randvec[2];
	float y = (not_cos_phi * u[0] * u[1]) * randvec[0] + (cos_phi + not_cos_phi * u[1] * u[1]) * randvec[1] + (-sin_phi * u[0]) * randvec[2];
	float z = (-sin_phi * u[1]) * randvec[0] + (sin_phi * u[0]) * randvec[1] + (cos_phi) * randvec[2];
	vec3 output = vec3(x, y, (b[2] / abs(b[2])) * z);
	return output;
}


vec3 relativeVector(vec3 macro_normal, vec3 offset) {
	vec3 a(0, 0, 1);
	vec3 b = macro_normal / macro_normal.length();
	//vec3 u = cross(a, b);
	vec3 u(-b[1], b[0], 0);
	float u_length = u.length();
	if (u_length == 0) { return a; }
	float sin_phi = min(u_length, 1.0f);
	u /= sin_phi;
	float phi = asin(sin_phi);
	float cos_phi = cos(phi);
	float not_cos_phi = 1 - cos_phi;

	float x = (cos_phi + not_cos_phi * u[0] * u[0]) * offset[0] + (not_cos_phi * u[0] * u[1]) * offset[1] + (sin_phi * u[1]) * offset[2];
	float y = (not_cos_phi * u[0] * u[1]) * offset[0] + (cos_phi + not_cos_phi * u[1] * u[1]) * offset[1] + (-sin_phi * u[0]) * offset[2];
	float z = (-sin_phi * u[1]) * offset[0] + (sin_phi * u[0]) * offset[1] + (cos_phi)*offset[2];
	/*vec3 output = vec3(x, y, (b[2] / abs(b[2])) * z);
	return output;*/
	if (b[2] < 0) {
		return vec3(x, y, -z);
	}
	return vec3(x, y, z);
}