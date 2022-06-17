#pragma once

#define sign(x) (x < 0 ? -1 : 1)

#include "hitable.h"

class grid : public hitable {
public:
	aabb grid_bbox;
	hitable** pointer_list;
	int num_hitables;
	int**** index_grid;
	int*** occupancy_grid;
	int* resolution;
	vec3 resolution_vec;

	grid() : resolution(0) {}
	grid(hitable** l, int n, float lambda = 4.0f) : pointer_list(l), num_hitables(n) {
		compute_bbox(l, n);
		vec3 d = grid_bbox._max - grid_bbox._min;
		float V = d[0] * d[1] * d[2];
		float cell_volume = lambda * n / V;
		float cube_root = cbrt(cell_volume);
		resolution_vec = d * cube_root;
		init(l, n);
	}
	grid(hitable** l, int n, int res) : pointer_list(l), num_hitables(n) {
		compute_bbox(l, n);
		resolution_vec = vec3(res, res, res);
		init(l, n);
	}
	grid(hitable** l, int n, vec3 res) : resolution_vec(res), pointer_list(l), num_hitables(n) {
		compute_bbox(l, n);
		init(l, n);
	}
	void init(hitable** l, int n) {
		resolution = new int[3];
		cout << "Grid resolution:  (";
		for (int i = 0; i < 3; i++) {
			resolution_vec[i] = max(1, int(resolution_vec[i]));
			resolution[i] = int(resolution_vec[i]);
			cout << resolution[i] << (i < 2 ? ", " : "");
		}
		cout << ")\n";

		occupancy_grid = new int** [resolution[0]];
		for (int i = 0; i < resolution[0]; i++) {
			occupancy_grid[i] = new int* [resolution[1]];
			for (int j = 0; j < resolution[1]; j++) {
				occupancy_grid[i][j] = new int[resolution[2]];
				for (int k = 0; k < resolution[2]; k++) {
					occupancy_grid[i][j][k] = 0;
				}
			}
		}

		aabb temp_box;
		int** min_points = new int* [n];
		int** max_points = new int* [n];
		for (int i = 0; i < n; i++) {
			min_points[i] = new int[3];
			max_points[i] = new int[3];
			if (l[i]->bounding_box(0, 0, temp_box)) {
				vec3 grid_diagonal = grid_bbox._max - grid_bbox._min;
				vec3 res_over_diagonal = resolution_vec * vec3(1, 1, 1) / grid_diagonal;
				vec3 rel_min = res_over_diagonal * (temp_box._min - grid_bbox._min);
				vec3 rel_max = res_over_diagonal * (temp_box._max - grid_bbox._min);
				for (int j = 0; j < 3; j++) {
					min_points[i][j] = max(int(rel_min[j]), 0);
					max_points[i][j] = min(int(rel_max[j]), resolution[j] - 1);
				}
				for (int x = min_points[i][0]; x <= max_points[i][0]; x++) {
					for (int y = min_points[i][1]; y <= max_points[i][1]; y++) {
						for (int z = min_points[i][2]; z <= max_points[i][2]; z++) {
							occupancy_grid[x][y][z]++;
						}
					}
				}
			}
			else {
				min_points[0][0] = -1;
			}
		}


		index_grid = new int*** [resolution[0]];
		for (int i = 0; i < resolution[0]; i++) {
			index_grid[i] = new int** [resolution[1]];
			for (int j = 0; j < resolution[1]; j++) {
				index_grid[i][j] = new int* [resolution[2]];
				for (int k = 0; k < resolution[2]; k++) {
					index_grid[i][j][k] = new int[occupancy_grid[i][j][k]];
					occupancy_grid[i][j][k] = 0;
				}
			}
		}


		for (int i = 0; i < n; i++) {
			if (min_points[0][0] != -1) {
				for (int x = min_points[i][0]; x <= max_points[i][0]; x++) {
					for (int y = min_points[i][1]; y <= max_points[i][1]; y++) {
						for (int z = min_points[i][2]; z <= max_points[i][2]; z++) {
							index_grid[x][y][z][occupancy_grid[x][y][z]++] = i;
						}
					}
				}
			}
		}
	}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const { box = grid_bbox; return true; };

	bool compute_bbox(hitable** list, int list_size);
};

bool grid::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
	hit_record temp_rec;
	bool hit_anything = false;
	double closest_so_far = t_max;

	if (grid_bbox.hit(r, t_min, t_max) == false) { return false; }

	float current_t = t_min;
	vec3 t0 = (grid_bbox._min - r.origin) / r.direction;
	vec3 t1 = (grid_bbox._max - r.origin) / r.direction;
	float lowest_t = min(t0[0], t1[0]);
	for (int i = 1; i < 3; i++) {
		float lower_t = min(t0[i], t1[i]);
		lowest_t = lower_t < lowest_t ? lower_t : lowest_t;
	}
	current_t = lowest_t > current_t ? lowest_t : current_t;


	vec3 cell_size = (grid_bbox._max - grid_bbox._min) / resolution_vec;
	vec3 direction_sign(sign(r.direction[0]), sign(r.direction[1]), sign(r.direction[2]));

	vec3 origin = r.point_at_parameter(current_t);
	vec3 rel_origin = (origin - grid_bbox._min) / cell_size;
	int cell[3];
	for (int i = 0; i < 3; i++) {
		cell[i] = max(min(int(rel_origin[i]), resolution[i] - 1), 0);
	}

	vec3 origin_cell(cell[0], cell[1], cell[2]);
	vec3 positive_direction = 0.5 * (direction_sign + vec3(1, 1, 1)); // because the next cell boundary in -ve direction is the start of the current cell

	vec3 next_t = (cell_size * (origin_cell + positive_direction) + grid_bbox._min - r.origin) / r.direction;
	vec3 delta_t = direction_sign * cell_size / r.direction;

	float* mailbox = new float[num_hitables](); // initialised to all zeros

	do {
		for (int i = 0; i < occupancy_grid[cell[0]][cell[1]][cell[2]]; i++) {
			int index = index_grid[cell[0]][cell[1]][cell[2]][i];
			if (mailbox[index] == 0.0f) { // if not tried intersection for that hitable yet
				mailbox[index] = 1.0f;
				if (pointer_list[index]->hit(r, t_min, closest_so_far, temp_rec)) {
					hit_anything = true;
					closest_so_far = temp_rec.t;
					rec = temp_rec;
				}
			}
		}

		int axis = 0;
		float min_next_t = next_t[0];
		if (next_t[1] < min_next_t) { axis = 1; min_next_t = next_t[1]; }
		if (next_t[2] < min_next_t) { axis = 2; min_next_t = next_t[2]; }

		current_t = min_next_t;
		cell[axis] += direction_sign[axis];
		next_t[axis] += delta_t[axis];
		vec3 temp = r.point_at_parameter(current_t);
		int i = 0;

	} while (
		(current_t < closest_so_far) 
		&& (min(cell[0], min(cell[1], cell[2])) >= 0) 
		&& (cell[0] < resolution[0]) && (cell[1] < resolution[1]) && (cell[2] < resolution[2])
		);

	delete[] mailbox;

	return hit_anything;
}


bool grid::compute_bbox(hitable** list, int list_size) {
	if (list_size < 1) {
		return false;
	}

	aabb temp_box;
	bool first_true = list[0]->bounding_box(0, 0, temp_box);
	if (!first_true) {
		return false;
	}
	else {
		grid_bbox = temp_box;
	}
	for (int i = 1; i < list_size; i++) {
		if (list[i]->bounding_box(0, 0, temp_box)) {
			grid_bbox = surrounding_box(grid_bbox, temp_box);
		}
		else {
			return false;
		}
	}
	return false;
}



//class octree : public hitable {
//public:
//	grid* root;
//
//	octree(hitable** l, int n, int objects_per_cell) {
//		root = new grid(l, n);
//	}
//
//	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
//	virtual bool bounding_box(float t0, float t1, aabb& box) const { box = grid_bbox; return true; };
//};