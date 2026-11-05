#include "student.h"

#include <algorithm>
#include <cmath>

// TASK 3a: a unit direction around n with density cos(theta) / pi.
vec3 cosine_sample_hemisphere(const vec3& n, float u1, float u2) {
  float r = std::sqrt(u1);
  float phi = 2 * PI * u2;
  vec3 t, b;
  make_basis(n, t, b);
  return normalize(t * (r * std::cos(phi)) + b * (r * std::sin(phi)) + n * std::sqrt(1 - u1));
}

// TASK 3b: the light reaching hit straight from one random point on the lamp.
vec3 sample_light(const Hit& hit, const PtScene& scene, Rng& rng) {
  const PtLight& light = scene.light;
  float u1 = rng.next(), u2 = rng.next();
  vec3 q = light.corner + light.e1 * u1 + light.e2 * u2;
  vec3 o = hit.p + hit.n * EPS;
  vec3 to_q = q - o;
  float d = length(to_q);
  vec3 l = to_q / d;
  float cos_x = dot(l, hit.n);
  float cos_l = -dot(l, light.normal());
  if(cos_x <= 0 || cos_l <= 0) {
    return vec3{0, 0, 0};
  }
  Hit blocker;
  if(pt_intersect(Ray{o, l}, scene, d - EPS, blocker)) {
    return vec3{0, 0, 0};
  }
  vec3 rho = pt_material(scene, hit).albedo;
  return rho / PI * light.emission * (cos_x * cos_l / (d * d) * light.area());
}

// TASK 3c: Russian roulette (EQ-10.08): false stops the path.
bool russian_roulette(vec3& throughput, int bounce, Rng& rng) {
  if(bounce < 3) {
    return true;
  }
  float q = std::min(0.95f, std::max(0.05f, 1 - luminance(throughput)));
  if(rng.next() < q) {
    return false;
  }
  throughput /= (1 - q);
  return true;
}

// TASK 3d: the light arriving along r, following one path.
vec3 path_trace(const Ray& r, const PtScene& scene, Rng& rng, const PtOptions& opt) {
  vec3 throughput{1, 1, 1}, result{0, 0, 0};
  Ray ray = r;
  for(int bounce = 0; bounce < opt.max_depth; ++bounce) {
    Hit hit;
    if(!pt_intersect(ray, scene, INFINITY, hit)) {
      break;
    }
    if(dot(hit.n, ray.d) > 0) {
      hit.n = -hit.n;
    }
    const PtMaterial& material = pt_material(scene, hit);
    if(luminance(material.emission) > 0) {
      if((bounce == 0 || !opt.nee) && dot(ray.d, scene.light.normal()) < 0) {
        result += throughput * material.emission;
      }
      break;
    }
    if(opt.nee) {
      result += throughput * sample_light(hit, scene, rng);
    }
    float u1 = rng.next(), u2 = rng.next();
    ray = Ray{hit.p + hit.n * EPS, cosine_sample_hemisphere(hit.n, u1, u2)};
    throughput *= material.albedo;
    if(opt.rr && !russian_roulette(throughput, bounce, rng)) {
      break;
    }
  }
  return result;
}
