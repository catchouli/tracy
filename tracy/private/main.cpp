#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <glm.hpp>
#include "stb_image_write.h"

const float PI = 3.1415926535897f;

struct Ray
{
  glm::vec3 o, d;
};

struct HitInfo
{
  glm::vec3 pos;
  glm::vec3 nrm;
};

struct Material
{
  enum Type { DIFFUSE, SPECULAR, DIELECTRIC };
  Type type;
  glm::vec3 albedo = glm::vec3(1.0f);
  glm::vec3 emissive = glm::vec3(0.0f);
};

class Object
{
public:
  Object(Material mat)
    : mat(mat) {}

  virtual bool intersect(const Ray& ray, float& outT) const = 0;
  virtual void getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const = 0;

  Material mat;
};

class Sphere : public Object
{
public:
  Sphere(glm::vec3 pos, float radius, Material mat)
    : pos(pos), rad(radius), Object(mat) {}

  bool intersect(const Ray& ray, float& outT) const override;
  void getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const override;

  glm::vec3 pos = glm::vec3(0.0f);
  float rad = 1.0f;
};

bool Sphere::intersect(const Ray& ray, float& outT) const
{
  glm::vec3 v = ray.o - pos;

  const float b = 2 * glm::dot(v, ray.d);
  const float c = glm::dot(v, v) - rad * rad;
  float delta = b * b - 4 * c;

  if (delta < 1e-4)
    return false;

  const float t1 = (-b - sqrt(delta)) / 2;
  const float t2 = (-b + sqrt(delta)) / 2;

  outT = (t1 < t2) ? t1 : t2; // get the first intersection only

  return true;
}

void Sphere::getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const
{
  hitInfo.pos = ray.o + ray.d * T;
  hitInfo.nrm = glm::normalize(hitInfo.pos - pos);
}

class Plane : public Object
{
public:
  Plane(glm::vec3 pos, glm::vec3 normal, Material mat)
   : pos(pos), nrm(normal), Object(mat) {}

  bool intersect(const Ray& ray, float& outT) const override;
  void getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const override;

  glm::vec3 pos;
  glm::vec3 nrm;
};

bool Plane::intersect(const Ray& ray, float& outT) const
{
  float denom = glm::dot(nrm, ray.d);
  if (abs(denom) > 0.0001f) // your favorite epsilon
  {
    float t = glm::dot((pos - ray.o), nrm) / denom;
    if (t >= 0)
    {
      outT = t;
      return true; // you might want to allow an epsilon here too
    }
  }
  return false;
}

void Plane::getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const
{
  hitInfo.pos = ray.o + ray.d * T;
  hitInfo.nrm = nrm;
}

Material diffuseWhite = { Material::DIFFUSE, glm::vec3(.75,.75,.75), glm::vec3() };
Material diffuseRed = { Material::DIFFUSE, glm::vec3(.75,.25,.25), glm::vec3() };
Material diffuseBlue = { Material::DIFFUSE, glm::vec3(.25,.25,.75), glm::vec3() };
Material diffuseGreen = { Material::DIFFUSE, glm::vec3(.25,.75,.25), glm::vec3() };
Material mirror = { Material::SPECULAR, glm::vec3(0.999f), glm::vec3() };
Material glass = { Material::DIELECTRIC, glm::vec3(0.999f), glm::vec3() };
Material light = { Material::DIFFUSE, glm::vec3(), glm::vec3(12,12,12) };

const Sphere spheres[] =
{
  Sphere { glm::vec3(-0.6f, -0.6f, 1.0f), 0.4f, glass },
  Sphere { glm::vec3(0.0f, -0.6f, 0.0f), 0.4f, mirror },
  Sphere { glm::vec3(0.6f, -0.6f, 1.0f), 0.4f, diffuseGreen },
  Sphere { glm::vec3(0.0f, 1.0 + 9.99f, 0.0f), 10.0f, light },
};

const Plane planes[] =
{
  Plane { glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, 0.0f, 1.0f), diffuseWhite }, // back
  Plane { glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), diffuseRed }, // left
  Plane { glm::vec3(1.0f,  0.0f, 0.0f), glm::vec3(-1.0f, 0.0f, 0.0f), diffuseBlue }, // right
  Plane { glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), diffuseWhite }, // bottom
  Plane { glm::vec3(0.0f,  1.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f), diffuseWhite }, // top
};

bool intersectScene(const Ray& ray, const Object*& hitObj, float& outT)
{
  hitObj = nullptr;
  outT = FLT_MAX;

  for (int i = 0; i < _countof(spheres); ++i)
  {
    const Sphere& sph = spheres[i];
    float T;
    if (sph.intersect(ray, T) && T > 0.0f && T < outT)
    {
      hitObj = &sph;
      outT = T;
    }
  }

  for (int i = 0; i < _countof(planes); ++i)
  {
    const Plane& plane = planes[i];
    float T;
    if (plane.intersect(ray, T) && T > 0.0f && T < outT)
    {
      hitObj = &plane;
      outT = T;
    }
  }

  return hitObj != nullptr;
}

float randf()
{
  return rand() / float(RAND_MAX);
}

// find perpendicular vectors for a normal
void createCoordinateSpace(const glm::vec3& n, glm::vec3& u, glm::vec3& v)
{
  // Plane equation
  //  A * x + B * y + C * z = D = 0
  // For the purpose of this, we aren't interested in D, and can set it to 0.
  // To find a perpendicular vector we can set y to 0, then we can use a cross product to get the other.
  // For normals where x < y we should use x = 0 instead.
  if (fabs(n.x) > fabs(n.y))
    u = glm::normalize(glm::vec3(n.x, 0.0f, -n.x));
  else
    u = glm::normalize(glm::vec3(0.0f, -n.z, n.y));
  v = glm::cross(n, u);
}

glm::vec3 uniformSampleHemisphere()
{
  const float r1 = randf();
  const float r2 = randf();
  float sinTheta = sqrtf(1.0f - r1 * r1);
  float phi = 2 * PI * r2;
  float x = sinTheta * cosf(phi);
  float z = sinTheta * sinf(phi);
  return glm::vec3(x, 1.0f, z);
}

glm::vec3 transformSample(const glm::vec3& sample, const glm::vec3& n, const glm::vec3& u, const glm::vec3& v)
{
  return glm::vec3(
    sample.x * v.x + sample.y * n.x + sample.z * u.x,
    sample.x * v.y + sample.y * n.y + sample.z * u.y,
    sample.x * v.z + sample.y * n.z + sample.z * u.z
  );
}

glm::vec3 radiance(const Ray& ray, int depth, int maxDepth)
{
  if (depth >= maxDepth)
    return glm::vec3();

  // Find nearest intersecting sphere
  const Object* nearest = nullptr;
  float nearestT = FLT_MAX;
  intersectScene(ray, nearest, nearestT);

  if (nearest == nullptr)
    return glm::vec3();

  glm::vec3 col = nearest->mat.albedo;

  // Russian roulette
  float p = glm::max(glm::max(nearest->mat.albedo.x, nearest->mat.albedo.y), nearest->mat.albedo.z);
  if (depth > 5)
  {
    if (randf() < p)
      col = col * (1.0f / p);
    else
      return nearest->mat.emissive;
  }

  HitInfo hitInfo;
  nearest->getHitInfo(ray, nearestT, hitInfo);

  if (nearest->mat.type == Material::DIFFUSE)
  {
    // Find coordinate space perpendicular to our normal
    glm::vec3 u, v;
    createCoordinateSpace(hitInfo.nrm, u, v);

    // Create a random sample about a hemisphere
    glm::vec3 sample = uniformSampleHemisphere();

    // Transform it into our normal's coordinate space
    sample = glm::normalize(transformSample(sample, hitInfo.nrm, u, v));

    Ray bounceRay = Ray { hitInfo.pos, sample };

    return nearest->mat.emissive + col * radiance(bounceRay, depth + 1, maxDepth);
  }
  else if (nearest->mat.type == Material::SPECULAR)
  {
    Ray reflectionRay = { hitInfo.pos, glm::reflect(ray.d, hitInfo.nrm) };
    return nearest->mat.emissive + col * radiance(reflectionRay, depth + 1, maxDepth);
  }
  else if (nearest->mat.type == Material::DIELECTRIC)
  {
    return glm::vec3();

    Ray reflectionRay = { hitInfo.pos, glm::reflect(ray.d, hitInfo.nrm) };
    glm::vec3 n1 = glm::dot(hitInfo.nrm, ray.d) < 0 ? hitInfo.nrm : -hitInfo.nrm;
    bool into = glm::dot(hitInfo.nrm, n1) > 0.0f;

    float nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = glm::dot(ray.d, n1), cos2t;

    // if total internal reflection, reflect
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // total internal reflection
      return nearest->mat.emissive + col * radiance(reflectionRay, depth + 1, maxDepth);

    // otherwise, choose reflection or refraction
    glm::vec3 tdir = glm::normalize(ray.d * nnt - hitInfo.nrm * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))));
    float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : glm::dot(tdir, hitInfo.nrm));
    float Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);

    Ray ray = Ray { ray.o, tdir };

    return nearest->mat.emissive + col * (depth > 2 ? (randf() < P ?   // Russian roulette 
      RP * radiance(reflectionRay, depth+1, maxDepth) : TP * radiance(ray, depth+1, maxDepth)) :
      Re * radiance(reflectionRay, depth+1, maxDepth) + Tr * radiance(ray, depth+1, maxDepth));
  }
  else
  {
    return glm::vec3();
  }
}

void trace(const Ray& ray, glm::vec4& rgba)
{
  const int spp = 40;
  const float avg = (1.0f / spp);

  glm::vec3 col;
  for (int i = 0; i < spp; ++i)
    col += avg * radiance(ray, 0, 10);

  rgba = glm::vec4(col, 1.0f);
}

int main(int argc, char** argv)
{
  const char* out = "out.png";
  const int width = 800;
  const int height = 600;
  const int comp = 4;

  const float aspect = width / (float)height; // assuming width > height 
  std::vector<uint8_t> pixels(width * height * comp);

  const float PI = 3.1415926535897f;

  float fov = 35.0f;

  Ray ray = { glm::vec3(0.0f, 0.0f, 5.0f) };

  float rowPercent = 0.0f;
  float rowPercentStep = 100.0f / height;
  float progressStep = 2;
  float lastProgressVal = 0;

  printf("progress: [");

  uint8_t* pix = pixels.data();
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      float px = (2.0f * ((x + 0.5f) / width) - 1.0f) * tan(fov / 2.0f * PI / 180.0f) * aspect;
      float py = (1.0f - 2.0f * ((y + 0.5f) / height)) * tan(fov / 2.0f * PI / 180.0f);

      ray.d = glm::normalize(glm::vec3(px, py, -1.0f));

      glm::vec4 rgba;

      trace(ray, rgba);

      rgba = glm::clamp(rgba, glm::vec4(0.0f), glm::vec4(1.0f));

      pix[0] = static_cast<uint8_t>(rgba.r * 255.0f);
      pix[1] = static_cast<uint8_t>(rgba.g * 255.0f);
      pix[2] = static_cast<uint8_t>(rgba.b * 255.0f);
      pix[3] = static_cast<uint8_t>(rgba.a * 255.0f);

      pix += comp;
    }

    rowPercent += rowPercentStep;
    if (rowPercent > lastProgressVal + progressStep)
    {
      printf("=");
      lastProgressVal += progressStep;
    }
  }

  printf("]\n");

  int res = stbi_write_png(out, width, height, comp, pixels.data(), width * comp * sizeof(uint8_t));

  if (res != 0) {
    printf("Successfully wrote image to %s\n", out);
  }
  else {
    printf("Failed to write image to %s\n", out);
  }
}