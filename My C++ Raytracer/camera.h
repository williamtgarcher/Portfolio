#ifndef CAMERAH
#define CAMERAH

#include "ray.h"
#include "random.h"

vec3 random_in_unit_disk() {
	vec3 p;
	do {
		p = 2.0 * vec3(random(), random(), 0) - vec3(1, 1, 0);
	} while (p.squared_length() >= 1.0);
	return p;
}

class camera {
public:
	vec3 u, v, w;
	float lens_radius;

	vec3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 origin;

	camera(vec3 position, vec3 view_direction, vec3 vup, float v_fov, float aspect, float aperture = 0.0, float focus_dist = 1.0, bool look_at = true) {
		float theta = v_fov * M_PI / 180;
		float half_height = tan(theta / 2);
		float half_width = aspect * half_height;
		lens_radius = aperture / 2;
		
		// If look_at is False, it uses view_direction as the vector the camera should point along. 
		// If look_at is True, it uses view_direction as a target coordinate to point the camera at.
		w = -unit_vector(look_at ? (view_direction - position) : view_direction);
		//w = unit_vector(position - view_direction);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);

		origin = position;
		lower_left_corner = origin - focus_dist * ((half_width * u) + (half_height * v) + w);
		horizontal = 2 * focus_dist * half_width * u;
		vertical = 2 * focus_dist * half_height * v;

		/*lower_left_corner = vec3(-half_width, -half_height, -1.0);
		horizontal = vec3(2 * half_width, 0.0, 0.0);
		vertical = vec3(0.0, 2 * half_height, 0.0);
		origin = vec3(0.0, 0.0, 0.0);*/
	}

	ray get_ray(float s, float t) {
		vec3 rd = lens_radius * random_in_unit_disk();
		//vec3 offset = rd.x() * u + rd.y() * v;
		vec3 offset = rd[0] * u + rd[1] * v;
		//vec3 offset_origin = origin + offset;
		//return ray(offset_origin, lower_left_corner + s * horizontal + t * vertical - offset_origin);
		return ray((origin + offset), lower_left_corner + s * horizontal + t * vertical - (origin + offset));
	}


	int* initial_media_inside;
	void init_media(hitable* scene) {
		initial_media_inside = new int[material::m];
		for (int i = 0; i < material::m; i++) { initial_media_inside[i] = 0; }
		hit_record temp_rec;
		temp_rec.t = 0;
		vec3 temp_vec;
		ray temp_ray;
		ray r_in = ray(origin, vec3(0, 0, 1)); // arbitrary direction
		r_in.media_inside = initial_media_inside;
		while (scene->hit(r_in, temp_rec.t + 0.0001, FLT_MAX, temp_rec)) {
			material::pointer_array[temp_rec.matptr_index]->scatter(r_in, temp_rec, temp_vec, temp_ray);
			r_in.media_inside = temp_ray.media_inside;
		}
		initial_media_inside = r_in.media_inside;
		for (int i = 0; i < material::m; i++) {
			initial_media_inside[i] = -initial_media_inside[i];
		}
	}

};

#endif