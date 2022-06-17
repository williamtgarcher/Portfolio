#pragma once

#include "texture.h"
#include "sphere.h"

class sky {
public:
	virtual vec3 value(const vec3& direction) const = 0;
};

class original_sky : public sky {
public:
	virtual vec3 value(const vec3& direction) const {
		vec3 unit_direction = unit_vector(direction);
		float t = 0.5 * (unit_direction.y() + 1.0);
		return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
	}
};

class constant_sky : public sky {
public:
	vec3 color;

	constant_sky() {}
	constant_sky(vec3 c) : color(c) {}
	virtual vec3 value(const vec3& direction) const { return color; }
};

class texture_sky : public sky {
public:
	texture* tex;

	texture_sky() {}
	texture_sky(texture* t) : tex(t) {}
	virtual vec3 value(const vec3& direction) const {
		float u, v;
		get_sphere_uv(unit_vector(direction), u, v);
		return tex->value(u, v, direction);
	}
};