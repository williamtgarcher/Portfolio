#define _USE_MATH_DEFINES

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <fstream>

#include <string>
#include <chrono>

#include "sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "material.h"
#include "random.h"
#include "perlin.h"
#include "constant_volume.h"
#include "material2.h"
#include "sky.h"
#include "grid.h"
#include "medium.h"

#include "random_scene.h"

using namespace std;


hitable* two_spheres() {
	hitable** list = new hitable * [2];
	texture* checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)), 10);
	list[0] = new sphere(vec3(0, -10, 0), 10, new lambertian(checker));
	list[1] = new sphere(vec3(0, 10, 0), 10, new lambertian(checker));
	return new hitable_list(list, 2);
}

hitable* perlin_spheres() {
	texture* pertext = new noise_texture(10.0);
	hitable** list = new hitable * [2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new sphere(vec3(0, 1, 0), 1, new lambertian(pertext));
	return new hitable_list(list, 2);
}

hitable* earth() {
	int px, py, pn;
	unsigned char* tex_data = stbi_load("2k_earth_daymap.jpg", &px, &py, &pn, 0);
	material* mat = new lambertian(new image_texture(tex_data, px, py));
	return new sphere(vec3(0, 0, 0), 2, mat);
}

hitable* simple_light() {
	hitable** list = new hitable * [3];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.4, 0.9, 0.9))));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(new constant_texture(vec3(0.9, 0.9, 0.4))));
	list[2] = new sphere(vec3(0, 7, 0), 2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
	return new hitable_list(list, 3);
}

hitable* simple_smoke() {
	hitable** list = new hitable * [3];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.4, 0.9, 0.9))));
	hitable* b = list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(new constant_texture(vec3(0.9, 0.9, 0.4))));
	list[1] = new constant_medium(b, 0.01, new constant_texture(vec3(0.5, 0.5, 0.5)));
	list[2] = new sphere(vec3(0, 7, 0), 2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
	return new hitable_list(list, 3);
}

hitable* material_test() {
	texture* pertext = new noise_texture(10.0);
	texture* black = new constant_texture(vec3(0, 0, 0));
	texture* pertext2 = new halves_texture(pertext, black, vec3(0, -1, 0));
	texture* white = new constant_texture(vec3(1, 1, 1));
	texture* red = new constant_texture(vec3(0.9, 0, 0));
	texture* flat = new constant_texture(vec3(0.5, 0.5, 1.0));
	texture* tiles = new checker_texture(new constant_texture(vec3(0.48, 0.48, 1.0)), new constant_texture(vec3(0.52, 0.52, 1.0)), 2);
	texture* mer = new constant_texture(vec3(0.0, 0.0, 0.05));
	hitable** list = new hitable * [2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new opaque(red, tiles, mer));
	list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
	return new hitable_list(list, 2);
}

hitable* transparent_material_test() {
	texture* pertext = new noise_texture(10.0);
	texture* black = new constant_texture(vec3(0, 0, 0));
	texture* white = new constant_texture(vec3(1, 1, 1));
	texture* red = new constant_texture(vec3(0.99, 0, 0));
	texture* flat = new constant_texture(vec3(0.5, 0.5, 1.0));
	texture* tiles = new checker_texture(new constant_texture(vec3(0.48, 0.48, 1.0)), new constant_texture(vec3(0.52, 0.52, 1.0)), 2);
	texture* mer = new constant_texture(vec3(0.7, 0.0, 0.05));
	hitable** list = new hitable * [2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new opaque(red, tiles, mer));
	list[1] = new sphere(vec3(0, 2, 0), 2, new translucent(white, flat, black, 1.5));
	return new hitable_list(list, 2);
}

hitable* grid_test() {
	hitable_list* scene = random_scene_list();
	hitable** list = new hitable * [2];
	int n = scene->list_size;
	list[0] = scene->list[0];
	hitable** scene_list = new hitable * [n - 1];
	for (int i = 0; i < n - 1; i++) { scene_list[i] = scene->list[i + 1]; }
	list[1] = new grid(scene_list, n - 1, 20);
	return new hitable_list(list, 2);
}

hitable* grid_test_dynamic() {
	hitable_list* scene = random_scene_list();
	return new grid(scene->list, scene->list_size);
}

hitable* medium_test() {
	hitable** list = new hitable * [3];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.5, 0.5, 0.5))));
	list[1] = new sphere(vec3(0, 2, 0), 10, new medium(vec3(1, 1, 1), 0.1, vec3(1, 1, 1), 0.01));
	list[2] = new sphere(vec3(-5, 2, 0), 2, new lambertian(new constant_texture(vec3(0.7, 0.7, 0.7))));
	return new hitable_list(list, 3);
}



int color_calls = 0; int total_color_calls = 0;
vec3 color(const ray& r, hitable* world, sky* skysphere, int depth) {
	color_calls++; total_color_calls++;

	hit_record rec;
	// t_min is slightly > 0 to prevent reflected rays colliding with the object they reflect from at very small t
	if (world->hit(r, 0.001, FLT_MAX, rec)) {
		ray scattered;
		vec3 attenuation;
		vec3 absorption = medium::media_absorption(r, rec);
		vec3 emitted = absorption * material::pointer_array[rec.matptr_index]->emitted(rec.u, rec.v, rec.p);
		emitted += medium::media_emission(r, rec);
		if (depth < 50 && material::pointer_array[rec.matptr_index]->scatter(r, rec, attenuation, scattered)) {
			if (rec.matptr_index >= material::m) { scattered.media_inside = r.media_inside; }
			return emitted + absorption * attenuation * color(scattered, world, skysphere, depth + 1);
		}
		else {
			return emitted;
		}
	}
	else {
		return skysphere->value(r.direction);
	}
}

int main() {
	chrono::system_clock::time_point tp = chrono::system_clock::now();
	chrono::system_clock::duration dtn = tp.time_since_epoch();
	string output_filename = "renders/Image [" + to_string(dtn.count()) + "].ppm";
	ofstream ImageFile("render.ppm");

	const int nx = 192*5, ny = 108*5;
	const int ns = 20;
	ImageFile << "P3\n" << nx << " " << ny << "\n255\n";

	
	sky* skysphere = new original_sky();
	//sky* skysphere = new constant_sky(vec3(0.1, 0.1, 0.1));

	

	hitable* world = medium_test();
	//hitable* world = grid_test_dynamic();
	//hitable* world = grid_test();
	//hitable* world = transparent_material_test();
	//hitable* world = material_test();
	//hitable* world = perlin_spheres();
	//hitable* world = random_scene();

	material::setup_pointer_array();

	srand(dtn.count());
	//srand(0);

	vec3 lookfrom(13, 2, 3);
	vec3 lookat(0, 1, 0);
	float dist_to_focus = 10.0;
	//float aperture = 0.1;
	float aperture = 0.0;

	camera cam(lookfrom, lookat, vec3(0, 1, 0), 20, float(nx) / float(ny), aperture, dist_to_focus);
	cam.init_media(world);

	int row_buffer[nx][3];
	int pixels_rendered = 0; int total_pixels_rendered = 0;
	for (int j = ny - 1; j >= 0; j--) {
		cerr << "\rScanlines remaining: " << j << " (" << to_string(float(100 * j) / float(ny)) << "%)" << ". Calls to color() per ray: " << to_string(float(color_calls) / float(ns * pixels_rendered)) << string(20, ' ') << flush;
		color_calls = 0; pixels_rendered = 0;
//#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nx; i++) {
			pixels_rendered++; total_pixels_rendered++;

			vec3 col(0, 0, 0);

			for (int s = 0; s < ns; s++) {
				float u = (float(i) + random()) / float(nx);
				float v = (float(j) + random()) / float(ny);

				ray r = cam.get_ray(u, v);
				r.media_inside = new int[material::m];
				for (int k = 0; k < material::m; k++) {
					r.media_inside[k] = cam.initial_media_inside[k];
				}

				col += color(r, world, skysphere, 0);
				delete[] r.media_inside;
			}

			col /= float(ns);

								
			row_buffer[i][0] = min(int(255.99 * sqrt(col[0])), 255);
			row_buffer[i][1] = min(int(255.99 * sqrt(col[1])), 255);
			row_buffer[i][2] = min(int(255.99 * sqrt(col[2])), 255);
		}
		for (int i = 0; i < nx; i++) {
			ImageFile << row_buffer[i][0] << " " << row_buffer[i][1] << " " << row_buffer[i][2] << "\n";
		}
	}

	
	ImageFile.close();
	cerr << "\rRender complete. Overall calls to color() per ray: " << to_string(float(total_color_calls) / float(ns * total_pixels_rendered)) << string(20, ' ') << flush;
	std::ifstream src("render.ppm", std::ios::binary);
	std::ofstream dst(output_filename, std::ios::binary);
	dst << src.rdbuf();
	return 0;
}