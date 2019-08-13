#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <glm.hpp>
#include "stb_image_write.h"
#include <random>

const float PI = 3.1415926535897f;
const float TWO_PI = (2.0f * PI);
const float INV_2PI = (1.0f / TWO_PI);

// Ray description
struct Ray
{
  glm::vec3 o, d;
};

// Type containing information about a ray hit
struct HitInfo
{
  glm::vec3 pos;
  glm::vec3 nrm;
};

// Material type
struct Material
{
  enum Type { DIFFUSE, SPECULAR, DIELECTRIC };
  Type type;
  glm::vec3 albedo = glm::vec3(1.0f);
  glm::vec3 emissive = glm::vec3(0.0f);
};

// Base object interface for intersection tests
class Object
{
public:
  Object(Material mat)
    : mat(mat) {}

  virtual bool intersect(const Ray& ray, float& outT) const = 0;
  virtual void getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const = 0;

  Material mat;
};

// A sphere object
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

// Sphere intersection - note I got this from somewhere but don't remember exactly where
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

// An infinite plane object
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

// Plane intersection
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool Plane::intersect(const Ray& ray, float& outT) const
{
  float denom = glm::dot(nrm, ray.d);
  if (abs(denom) > 1e-6)
  {
    float t = glm::dot((pos - ray.o), nrm) / denom;
    if (t >= 1e-6) {
      outT = t;
      return true;
    }
  }
  return false;
}

void Plane::getHitInfo(const Ray& ray, const float T, HitInfo& hitInfo) const
{
  hitInfo.pos = ray.o + ray.d * T;
  hitInfo.nrm = nrm;
}

// Disk object
class Disk : public Plane
{
public:
  Disk(glm::vec3 pos, glm::vec3 normal, float radius, Material mat)
   : rad(radius), Plane(pos, normal, mat) {}

  bool intersect(const Ray& ray, float& outT) const override;

  float rad;
};

// Disk intersection
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool Disk::intersect(const Ray& ray, float& outT) const
{
  float t = 0.0f;
  if (!Plane::intersect(ray, t))
    return false;

  glm::vec3 p = ray.o + ray.d * t;
  glm::vec3 v = p - pos;
  float d2 = glm::dot(v, v);

  if (d2 <= rad * rad)
  {
    outT = t;
    return true;
  }
  else
  {
    return false;
  }
}


// Scene definition
struct Scene
{
  std::vector<std::shared_ptr<Object>> objects;
};

// Generates a random float from 0..1
float randf()
{
  static std::default_random_engine e;
  static std::uniform_real_distribution<float> dis(0, 1); // rage 0 - 1
  return dis(e);
}

// Tests intersections against a scene description and returns
// the hit object and distance down the ray as hitObj and outT
void intersectScene(const Ray& ray, const Scene& scene, const Object*& hitObj, float& outT)
{
  hitObj = nullptr;
  outT = FLT_MAX;

  for (auto it = scene.objects.begin(); it != scene.objects.end(); ++it)
  {
    const Object& obj = **it;
    float T;
    if (obj.intersect(ray, T) && T > 0.0f && T < outT)
    {
      hitObj = &obj;
      outT = T;
    }
  }
}

// Find perpendicular vectors for a normal to construct a coordinate space about it
void createCoordinateSpace(const glm::vec3& n, glm::vec3& u, glm::vec3& v)
{
  // Plane equation
  //  A * x + B * y + C * z = D = 0
  // For the purpose of this, we aren't interested in D, and can set it to 0.
  // To find a perpendicular vector we can set y to 0, then we can use a cross product to get the other.
  // For normals where x < y we should use x = 0 instead.
  if (fabs(n.x) > fabs(n.y))
    u = glm::vec3(n.z, 0.0f, -n.x) / sqrtf(n.x * n.x + n.z * n.z);
  else
    u = glm::vec3(0.0f, -n.z, n.y) / sqrtf(n.y * n.y + n.z * n.z);
  v = glm::cross(n, u);
}

// Transforms a vector into a coordinate space obtained from createCoordinateSpace
glm::vec3 transformIntoSpace(const glm::vec3& sample, const glm::vec3& n, const glm::vec3& u, const glm::vec3& v)
{
  return glm::vec3(
    sample.x * v.x + sample.y * n.x + sample.z * u.x,
    sample.x * v.y + sample.y * n.y + sample.z * u.y,
    sample.x * v.z + sample.y * n.z + sample.z * u.z
  );
}

// Generate a uniform hemisphere sample
void uniformHemisphereSample(glm::vec3& dir, float& pdf)
{
  float u1 = randf(), u2 = randf();

  const float r = sqrt(u1);
  const float theta = 2 * PI * u2;

  const float x = r * cos(theta);
  const float y = r * sin(theta);

  dir = glm::vec3(x, y, sqrt(glm::max(0.0f, 1 - u1)));
  pdf = INV_2PI;

  //// Random (uniformly distributed) point on a hemisphere based on two uniform random floats 0..1
  //// http://www.rorydriscoll.com/2009/01/07/better-sampling/
  //// Not sure if this actually produces a uniform distribution but it looks better than the function I originally wrote
  //float u1 = randf(), u2 = randf();
  //const float r = sqrt(1.0f - u1 * u1);
  //const float phi = TWO_PI * u2;
  //dir = glm::vec3(cos(phi) * r, sin(phi) * r, u1);

  //// The probability density function for a uniform distribution is constant
  //pdf = INV_2PI;
}

// Calculates the ratio for a lambertian brdf
// https://computergraphics.stackexchange.com/a/4280
float lambertBRDF(float cosTheta, const glm::vec3& p)
{
  return cosTheta / PI;
}

// Calculate incoming radiance along ray
glm::vec3 radiance(const Ray& ray, const Scene& scene, int depth, int maxDepth)
{
  if (depth >= maxDepth)
    return glm::vec3();

  // Find nearest scene intersection
  const Object* nearest = nullptr;
  float nearestT = FLT_MAX;
  intersectScene(ray, scene, nearest, nearestT);

  // Just return black if we didn't hit anything
  if (nearest == nullptr)
    return glm::vec3();

  // The color from the material of the hit object
  glm::vec3 col = nearest->mat.albedo;

  // Russian roulette - at a certain depth we want to start randomly discarding samples based on the current radiance
  float maxRadiance = glm::max(glm::max(nearest->mat.albedo.x, nearest->mat.albedo.y), nearest->mat.albedo.z);
  if (depth > 5)
  {
    if (randf() < maxRadiance)
      col = col * (1.0f / maxRadiance);
    else
      return nearest->mat.emissive;
  }

  // Get the hit info for this hit
  HitInfo hitInfo;
  nearest->getHitInfo(ray, nearestT, hitInfo);

  // Compute lambertian materials
  // https://computergraphics.stackexchange.com/a/4280
  if (nearest->mat.type == Material::DIFFUSE)
  {
    // Correct normal for ray direction
    glm::vec3 n = glm::dot(hitInfo.nrm, ray.d) < 0 ? hitInfo.nrm : -hitInfo.nrm;

    // Find coordinate space perpendicular to our normal
    glm::vec3 u, v;
    createCoordinateSpace(n, u, v);

    // Create a random sample about a hemisphere
    glm::vec3 sampleDir;
    float pdf;
    uniformHemisphereSample(sampleDir, pdf);

    // Transform sample dir from tangent space to world
    sampleDir = sampleDir.z * n + sampleDir.x * u + sampleDir.y * v;

    // Compute ncoming radiance
    glm::vec3 incoming = radiance(Ray{ hitInfo.pos, sampleDir }, scene, depth + 1, maxDepth);

    // Compute brdf
    // https://computergraphics.stackexchange.com/questions/4277/derivation-of-wikipedias-path-tracing-diffuse-brdf/4280#4280
    glm::vec3 outgoing = (col / PI) * incoming * glm::dot(sampleDir, n) / INV_2PI;

    // Add emissive term
    // TODO: cosine distribution and / pdf?
    return nearest->mat.emissive + outgoing;
  }
  // Compute specular materials
  else if (nearest->mat.type == Material::SPECULAR)
  {
    Ray reflectionRay = { hitInfo.pos, glm::reflect(ray.d, hitInfo.nrm) };
    return nearest->mat.emissive + col * radiance(reflectionRay, scene, depth + 1, maxDepth);
  }
  // Compute dialectric materials (currently a disaster area)
  else if (nearest->mat.type == Material::DIELECTRIC)
  {
    return glm::vec3();

    const float airIOR = 1.0f;
    const float glassIOR = 1.5f;

    // whether we're entering or exiting the glass
    bool entering = glm::dot(hitInfo.nrm, ray.d) < 0;

    // the corrected normal based on which direction we're going through the medium
    glm::vec3 n = entering ? hitInfo.nrm : -hitInfo.nrm;

    // the ratio of indices of reflection based on direction
    float eta = entering ? (airIOR / glassIOR) : (glassIOR / airIOR);

    // dot product of the ray direction and the normal
    float rayCos = glm::dot(ray.d, n);

    // Total internal reflection
    // TODO: no idea if this works, the code never seems to get hit
    if (eta < 1.0f)
    {
      float rayAngle = acos(rayCos);
      float criticalAngle = asin(eta);
      if (rayAngle > criticalAngle)
      {
        Ray reflectionRay = { hitInfo.pos, glm::reflect(ray.d, hitInfo.nrm) };
        return nearest->mat.emissive + col * radiance(reflectionRay, scene, depth + 1, maxDepth);
      }
    }

    // Give up now, we need to handle remaining reflection/refraction still
    return glm::vec3();
  }
  else
  {
    return glm::vec3();
  }
}

// Trace a certain number of samples along a ray in a given scene and return the radiance as rgba
glm::vec4 trace(const Ray& ray, const Scene& scene, int spp)
{
  const int maxDepth = 10;
  const float avg = (1.0f / spp);

  glm::vec3 col;
  for (int i = 0; i < spp; ++i)
    col += avg * radiance(ray, scene, 0, maxDepth);

  return glm::vec4(col, 1.0f);
}

// Materials
//                        Material type         Albedo                          Emissive
Material diffuseWhite = { Material::DIFFUSE,    glm::vec3(0.75f),               glm::vec3() };
Material diffuseRed =   { Material::DIFFUSE,    glm::vec3(0.75f, 0.25f, 0.25f), glm::vec3() };
Material diffuseRedW =  { Material::DIFFUSE,    glm::vec3(0.75f, 0.25f, 0.25f), glm::vec3(0.18f) };
Material diffuseBlue =  { Material::DIFFUSE,    glm::vec3(0.25f, 0.25f, 0.75f), glm::vec3() };
Material diffuseGreen = { Material::DIFFUSE,    glm::vec3(0.25f, 0.75f, 0.25f), glm::vec3() };
Material mirror =       { Material::SPECULAR,   glm::vec3(0.999f),              glm::vec3() };
Material glass =        { Material::DIELECTRIC, glm::vec3(0.999f),              glm::vec3() };
Material light =        { Material::DIFFUSE,    glm::vec3(),                    glm::vec3(12.0f) };

Scene cornellBox =
{
  {
    std::make_shared<Disk>(glm::vec3(0.0f, 1.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f), 0.25f, light),

    std::make_shared<Sphere>(glm::vec3(-0.6f, -0.6f, 1.0f), 0.4f, glass),
    std::make_shared<Sphere>(glm::vec3(0.0f, -0.6f, 0.0f), 0.4f, mirror),
    std::make_shared<Sphere>(glm::vec3(0.6f, -0.6f, 1.0f), 0.4f, diffuseGreen),
    //std::make_shared<Sphere>(glm::vec3(0.0f, 0.0f, 0.0f), 0.5f, light),
    //std::make_shared<Sphere>(glm::vec3(0.0f, 1.0f + 9.99f, 0.0f), 10.0f, light),

    std::make_shared<Plane>(glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, 0.0f, 1.0f), diffuseWhite), // back
    //std::make_shared<Plane>(glm::vec3(0.0f, 0.0f, 1.0f),  glm::vec3(0.0f, 0.0f, -1.0f), diffuseWhite), // front
    std::make_shared<Plane>(glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), diffuseRed), // left
    std::make_shared<Plane>(glm::vec3(1.0f,  0.0f, 0.0f), glm::vec3(-1.0f, 0.0f, 0.0f), diffuseBlue), // right
    std::make_shared<Plane>(glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), diffuseWhite), // bottom
    std::make_shared<Plane>(glm::vec3(0.0f,  1.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f), diffuseWhite), // top
  },
};

Scene leftWall =
{
  {
    std::make_shared<Plane>(glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), diffuseRedW), // left wall
  }
};

Scene furnace =
{
  {
    //std::make_shared<Sphere>(glm::vec3(0.0f), 1.0f, Material{ Material::DIFFUSE, glm::vec3(1.0f), glm::vec3() }),
    std::make_shared<Sphere>(glm::vec3(0.0f), 1000.0f, Material{ Material::DIFFUSE, glm::vec3(), glm::vec3(0.18f) }),

    std::make_shared<Plane>(glm::vec3(-2.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), Material{ Material::DIFFUSE, glm::vec3(), glm::vec3(0.18f) }), // top
  }
};

// Entry point - creates a pixel buffer and passes rays to trace()
int main(int argc, char** argv)
{
  const char* out = "out.png";
  const int width = 800;
  const int height = 600;
  const int comp = 4;

  const float aspect = width / (float)height; // assuming width > height 
  std::vector<uint8_t> pixels(width * height * comp);

  // Render parameters
  const Scene& sceneToRender = cornellBox; // The scene to render
  const glm::vec3 camPos = glm::vec3(0.0f, 0.0f, 5.0f); // The camera position
  const float fov = 35.0f; // The fov in degrees
  const int spp = 500; // The number of samples to take per pixel
  const int supersamples = 1; // Extra samples to take in each axis

  // Housekeeping for the progress meter
  float rowPercent = 0.0f;
  float rowPercentStep = 100.0f / height;
  float progressStep = 2;
  float lastProgressVal = 0;

  // Print first part of progress display
  printf("progress: [");

  // Coefficient for the fov
  float fovCoeff = tan(fov / 2.0f * PI / 180.0f);

  // Iterate over each pixel and trace scene
  uint8_t* pix = pixels.data();
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      glm::vec4 rgba = glm::vec4();

      for (int sx = 0; sx < supersamples; ++sx)
      {
        for (int sy = 0; sy < supersamples; ++sy)
        {
          const float ss = 1.0f / supersamples;
          const float ss2 = 1.0f / (supersamples * supersamples);

          // Construct ray based on fov and pixel
          float px = (2.0f * ((x + sx * ss) / width) - 1.0f) * fovCoeff * aspect;
          float py = (1.0f - 2.0f * ((y + sy * ss) / height)) * fovCoeff;
          Ray ray = Ray{ camPos, glm::normalize(glm::vec3(px, py, -1.0f)) };

          // Trace ray
          rgba += ss2 * trace(ray, sceneToRender, spp);
        }
      }

      // Convert to int
      pix[0] = static_cast<uint8_t>(rgba.r * 255.0f);
      pix[1] = static_cast<uint8_t>(rgba.g * 255.0f);
      pix[2] = static_cast<uint8_t>(rgba.b * 255.0f);
      pix[3] = static_cast<uint8_t>(rgba.a * 255.0f);

      pix += comp;
    }

    // Update progress
    rowPercent += rowPercentStep;
    if (rowPercent > lastProgressVal + progressStep)
    {
      printf("=");
      lastProgressVal += progressStep;
    }
  }

  // Finish progress
  printf("]\n");

  // Write image
  int res = stbi_write_png(out, width, height, comp, pixels.data(), width * comp * sizeof(uint8_t));

  if (res != 0) {
    printf("Successfully wrote image to %s\n", out);
  }
  else {
    printf("Failed to write image to %s\n", out);
  }
}