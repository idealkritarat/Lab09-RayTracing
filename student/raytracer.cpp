#include "student.h"

#include <algorithm>
#include <cmath>

// TASK 1a: the ray from the camera through the centre of pixel (x, y).
Ray camera_ray(int x, int y, int width, int height, float fov) {
  float scale = std::tan(fov / 2);
  float sx = (2 * (float(x) + 0.5f) / float(width) - 1) * scale * float(width) / float(height);
  float sy = -(2 * (float(y) + 0.5f) / float(height) - 1) * scale;
  return Ray{vec3{0, 0, 0}, normalize(vec3{sx, sy, -1})};
}

// TASK 1b: the smallest t > EPS where r meets the sphere (EQ-10.02).
bool ray_sphere(const Ray& r, const Sphere& s, float& t) {
  vec3 oc = r.o - s.center;
  float a = dot(r.d, r.d);
  float b = 2 * dot(r.d, oc);
  float c = dot(oc, oc) - s.radius * s.radius;
  float disc = b * b - 4 * a * c;
  if(disc < 0) {
    return false;
  }
  float root = std::sqrt(disc);
  float near_t = (-b - root) / (2 * a);
  float far_t = (-b + root) / (2 * a);
  if(near_t > EPS) {
    t = near_t;
    return true;
  }
  if(far_t > EPS) {
    t = far_t;
    return true;
  }
  return false;
}

// TASK 1c: d reflected about the unit normal n (EQ-10.04).
vec3 reflect(const vec3& d, const vec3& n) {
  return d - n * (2 * dot(d, n));
}

// TASK 1c: d refracted through a surface whose normal n points out of the object.
// false on total internal reflection.
bool refract(const vec3& d, const vec3& n, float ior, vec3& out) {
  float cos_i = -dot(d, n);
  float eta = 1 / ior;
  vec3 facing = n;
  if(cos_i < 0) {
    cos_i = -cos_i;
    eta = ior;
    facing = -n;
  }
  float k = 1 - eta * eta * (1 - cos_i * cos_i);
  if(k < 0) {
    return false;
  }
  out = d * eta + facing * (eta * cos_i - std::sqrt(k));
  return true;
}

static Ray bounce_ray(const Hit& hit, const vec3& dir) {
  vec3 offset = hit.n * (dot(dir, hit.n) < 0 ? -EPS : EPS);
  return Ray{hit.p + offset, dir};
}

// TASK 1d: the colour seen along r: Phong, shadows, reflection, refraction (Whitted).
vec3 cast_ray(const Ray& r, const Scene& scene, int depth) {
  Hit hit;
  if(depth > MAX_DEPTH || !scene_intersect(r, scene, hit)) {
    return scene.background;
  }
  const Material& m = *hit.material;
  if(STAGE < 2) {
    return m.diffuse;
  }

  vec3 reflected{0, 0, 0}, refracted{0, 0, 0};
  if(STAGE >= 5) {
    reflected = cast_ray(bounce_ray(hit, normalize(reflect(r.d, hit.n))), scene, depth + 1);
  }
  vec3 refract_dir;
  if(STAGE >= 6 && refract(r.d, hit.n, m.ior, refract_dir)) {
    refracted = cast_ray(bounce_ray(hit, normalize(refract_dir)), scene, depth + 1);
  }

  float diffuse = 0, specular = 0;
  for(const Light& light : scene.lights) {
    vec3 to_light = light.position - hit.p;
    float dist = length(to_light);
    vec3 l = to_light / dist;
    if(STAGE >= 4) {
      Hit blocker;
      if(scene_intersect(bounce_ray(hit, l), scene, blocker) && blocker.t < dist) {
        continue;
      }
    }
    diffuse += light.intensity * std::max(0.0f, dot(l, hit.n));
    if(STAGE >= 3) {
      specular += light.intensity * std::pow(std::max(0.0f, -dot(reflect(-l, hit.n), r.d)), m.specular_exponent);
    }
  }

  return m.diffuse * (diffuse * m.albedo[0]) + vec3{1, 1, 1} * (specular * m.albedo[1]) +
         reflected * m.albedo[2] + refracted * m.albedo[3];
}
