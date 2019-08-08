#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <glm.hpp>
#include "stb_image_write.h"
#include "math.h"

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
  enum Type { DIFFUSE, SPECULAR };
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
Material mirror = { Material::SPECULAR, glm::vec3(0.99f), glm::vec3() };
Material light = { Material::DIFFUSE, glm::vec3(), glm::vec3(12,12,12) };

const Sphere spheres[] =
{
  Sphere { glm::vec3(0.4f, -0.6f, 1.0f), 0.4f, diffuseGreen },
  Sphere { glm::vec3(-0.4f, -0.6f, 0.0f), 0.4f, mirror },
  //Sphere { glm::vec3(0.0f, 1.0f, 0.0f), 0.4f, light },
  //Sphere { glm::vec3(50.0f, 40.8f, 81.6f),   44.2,  glm::vec3(0.0f), glm::vec3(1.0f, 0.0f, 1.0f) },
  //Sphere { glm::vec3(1e5 + 1, 40.8, 81.6),   1e5,   glm::vec3(0.0f), glm::vec3(.75,.25,.25)      }, // left
  //Sphere { glm::vec3(-1e5 + 99, 40.8, 81.6), 1e5,   glm::vec3(0.0f), glm::vec3(.25,.25,.75)      }, // right
  //Sphere { glm::vec3(50, 40.8, 1e5),         1e5f,   glm::vec3(0.0f), glm::vec3(.75,.75,.75)      }, // back
  //Sphere { glm::vec3(50, 40.8, -1e5 + 170),  1e5f,   glm::vec3(0.0f), glm::vec3(0.0f)             }, // front
  //Sphere { glm::vec3(50, 1e5, 81.6),         1e5f,   glm::vec3(0.0f), glm::vec3(.75,.75,.75)      }, // bottom
  //Sphere { glm::vec3(50, -1e5 + 81.6, 81.6), 1e5f,   glm::vec3(0.0f), glm::vec3(.75,.75,.75)      }, // top
  //Sphere { glm::vec3(27, 16.5, 47),          16.5f,  glm::vec3(0.0f), 0.999f * glm::vec3(1,1,1)   }, // mirror
  //Sphere { glm::vec3(73, 16.5, 78),          16.5f,  glm::vec3(0.0f), 0.999f * glm::vec3(1,1,1)   }, // glass
  //Sphere { glm::vec3(50, 681.6-0.27, 81.6),  600.0f, glm::vec3(12.0f, 12.0f, 12.0f), glm::vec3(0.0f), } // light
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

glm::vec3 traceImpl(const Ray& ray, int depth, int maxDepth)
{
  glm::vec3 col;

  // A directional light
  const glm::vec3 lightPos = glm::vec3(0.0f, 0.7f, 0.0f);

  // Find nearest intersecting sphere
  const Object* nearest = nullptr;
  float nearestT = FLT_MAX;
  intersectScene(ray, nearest, nearestT);

  if (nearest != nullptr)
  {
    HitInfo hitInfo;
    nearest->getHitInfo(ray, nearestT, hitInfo);

    glm::vec3 L = glm::normalize(lightPos - hitInfo.pos);

    // Shadow ray
    Ray shadowRay = { hitInfo.pos, L };
    float lightDist = glm::length(hitInfo.pos - lightPos);
    const Object* shadowHit; float shadowT;
    intersectScene(shadowRay, shadowHit, shadowT);
    bool inShadow = shadowHit != nullptr && shadowT > 0.0f && shadowT < lightDist;

    // point light
    float ratio = inShadow ? 0.0f : glm::clamp(glm::dot(hitInfo.nrm, L), 0.0f, 1.0f);

    col += glm::vec3(nearest->mat.emissive + ratio * nearest->mat.albedo);

    if (nearest->mat.type == Material::SPECULAR && depth < maxDepth)
    {
      Ray reflectionRay = { hitInfo.pos, glm::reflect(ray.d, hitInfo.nrm) };
      col += 0.5f * traceImpl(reflectionRay, depth+1, maxDepth);
    }
  }

  return col;
}

void trace(const Ray& ray, glm::vec4& rgba)
{
  glm::vec3 col = traceImpl(ray, 0, 5);
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
  }

  int res = stbi_write_png(out, width, height, comp, pixels.data(), width * comp * sizeof(uint8_t));

  if (res != 0) {
    printf("Successfully wrote image to %s\n", out);
  }
  else {
    printf("Failed to write image to %s\n", out);
  }
}