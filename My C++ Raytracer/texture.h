#pragma once

#include "ray.h"

class texture {
public:
	virtual vec3 value(float u, float v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
	vec3 color;

	constant_texture() {}
	constant_texture(vec3 c) : color(c) {}
	virtual vec3 value(float u, float v, const vec3& p) const {
		return color;
	}
};

class checker_texture : public texture {
public:
	texture* odd;
	texture* even;
	const float f = 10.0;

	checker_texture() {}
	checker_texture(texture* t0, texture* t1) : even(t0), odd(t1) {}
	checker_texture(texture* t0, texture* t1, float freq): even(t0), odd(t1), f(freq) {}
	virtual vec3 value(float u, float v, const vec3& p) const {
		float sines = sin(f * p.x()) * sin(f * p.y()) * sin(f * p.z());
		if (sines < 0) {
			return odd->value(u, v, p);
		}
		else {
			return even->value(u, v, p);
		}
	}
};

class halves_texture : public texture {
public:
	texture* top;
	texture* bottom;
	vec3 axis;

	halves_texture() {}
	halves_texture(texture* t0, texture* t1, vec3 v) : top(t0), bottom(t1), axis(unit_vector(v)) {}
	virtual vec3 value(float u, float v, const vec3& p) const {
		float phi = 2 * M_PI * (1 - u) - M_PI;
		float theta = M_PI * v - M_PI / 2;
		float x = cos(theta) * cos(phi);
		float z = cos(theta) * sin(phi);
		float y = sin(theta);

		if (dot(vec3(x,y,z), axis) < 0) {
			return bottom->value(u, v, p);
		}
		else {
			return top->value(u, v, p);
		}
	}
};

class image_texture : public texture {
public:
	unsigned char* data;
	int nx, ny;

	image_texture() {}
	image_texture(unsigned char* pixels, int A, int B) : data(pixels), nx(A), ny(B) {}
	virtual vec3 value(float u, float v, const vec3& p) const;
};

vec3 image_texture::value(float u, float v, const vec3& p) const {
	int i = u * nx;
	int j = (1 - v) * ny - 0.001;

	i = (i < 0) ? 0 : i;
	j = (j < 0) ? 0 : j;

	i = (i > nx - 1) ? nx - 1 : i;
	j = (j > nx - 1) ? nx - 1 : j;

	int index = 3 * i + 3 * nx * j;
	float r = int(data[index]) / 255.0f;
	float g = int(data[index + 1]) / 255.0f;
	float b = int(data[index + 2]) / 255.0f;

	return vec3(r, g, b);
}