#include "student.h"

#include <algorithm>
#include <cmath>

// TASK 2a: Moller-Trumbore (DEF-10.01): does r hit the triangle at t > EPS? Write t, u, v.
bool ray_triangle(const Ray& r, const vec3& v0, const vec3& v1, const vec3& v2,
                  float& t, float& u, float& v) {
  vec3 e1 = v1 - v0, e2 = v2 - v0;
  vec3 h = cross(r.d, e2);
  float a = dot(e1, h);
  if(std::fabs(a) < 1e-8f) {
    return false;
  }
  vec3 s = r.o - v0;
  u = dot(s, h) / a;
  if(u < 0 || u > 1) {
    return false;
  }
  vec3 q = cross(s, e1);
  v = dot(r.d, q) / a;
  if(v < 0 || u + v > 1) {
    return false;
  }
  t = dot(e2, q) / a;
  return t > EPS;
}

static void test_triangle(const Ray& r, const Mesh& m, int index, float& best, Hit& hit, bool& found) {
  const auto& tri = m.tris[size_t(index)];
  const vec3& v0 = m.positions[size_t(tri[0])];
  const vec3& v1 = m.positions[size_t(tri[1])];
  const vec3& v2 = m.positions[size_t(tri[2])];
  float t = 0, u = 0, v = 0;
  count_triangle_test();
  if(!ray_triangle(r, v0, v1, v2, t, u, v) || t >= best) {
    return;
  }
  best = t;
  found = true;
  hit.t = t;
  hit.tri = index;
  if(!m.normals.empty()) {
    hit.n = normalize(m.normals[size_t(tri[0])] * (1 - u - v) + m.normals[size_t(tri[1])] * u +
                      m.normals[size_t(tri[2])] * v);
  } else {
    hit.n = normalize(cross(v1 - v0, v2 - v0));
  }
}

// TASK 2b: the closest triangle of m below t_max, testing every one (count_triangle_test()).
// Fill hit.t, hit.tri and hit.n.
bool intersect_mesh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  float best = t_max;
  bool found = false;
  for(int i = 0; i < int(m.tris.size()); ++i) {
    test_triangle(r, m, i, best, hit, found);
  }
  return found;
}

// TASK 2c (optional): the slab test (EQ-10.03). inv_d is 1 / r.d; write the entry distance.
bool ray_aabb(const Ray& r, const vec3& inv_d, const AABB& box, float t_max, float& t_entry) {
  float entry = -INFINITY, exit = INFINITY;
  for(int k = 0; k < 3; ++k) {
    float t0 = (box.lo[k] - r.o[k]) * inv_d[k];
    float t1 = (box.hi[k] - r.o[k]) * inv_d[k];
    entry = std::max(entry, std::min(t0, t1));
    exit = std::min(exit, std::max(t0, t1));
  }
  t_entry = entry;
  return entry <= exit && exit > 0 && entry < t_max;
}

// TASK 2d (optional): TASK 2b through m.bvh (count_node_visit() for every node).
bool intersect_bvh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  if(m.bvh.nodes.empty()) {
    return false;
  }
  vec3 inv_d{1 / r.d.x, 1 / r.d.y, 1 / r.d.z};
  float best = t_max;
  bool found = false;
  int stack[128];
  int top = 0;
  stack[top++] = 0;
  while(top > 0) {
    int index = stack[--top];
    const BVHNode& node = m.bvh.nodes[size_t(index)];
    count_node_visit();
    float entry = 0;
    if(!ray_aabb(r, inv_d, node.box, best, entry)) {
      continue;
    }
    if(node.count > 0) {
      for(int i = node.first; i < node.first + node.count; ++i) {
        test_triangle(r, m, m.bvh.order[size_t(i)], best, hit, found);
      }
      continue;
    }
    int near_child = index + 1, far_child = node.right;
    float near_t = INFINITY, far_t = INFINITY;
    if(!ray_aabb(r, inv_d, m.bvh.nodes[size_t(near_child)].box, best, near_t)) {
      near_t = INFINITY;
    }
    if(!ray_aabb(r, inv_d, m.bvh.nodes[size_t(far_child)].box, best, far_t)) {
      far_t = INFINITY;
    }
    if(far_t < near_t) {
      std::swap(near_child, far_child);
    }
    stack[top++] = far_child;
    stack[top++] = near_child;
  }
  return found;
}
