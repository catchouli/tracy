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

// Plane intersection - note I got this from somewhere but don't remember exactly where
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

// Scene definition
struct Scene
{
  std::vector<Sphere> spheres;
  std::vector<Plane> planes;
};

// Materials
//                        Material type         Albedo                          Emissive
Material diffuseWhite = { Material::DIFFUSE,    glm::vec3(0.75f),               glm::vec3() };
Material diffuseRed =   { Material::DIFFUSE,    glm::vec3(0.75f, 0.25f, 0.25f), glm::vec3() };
Material diffuseBlue =  { Material::DIFFUSE,    glm::vec3(0.25f, 0.25f, 0.75f), glm::vec3() };
Material diffuseGreen = { Material::DIFFUSE,    glm::vec3(0.25f, 0.75f, 0.25f), glm::vec3() };
Material mirror =       { Material::SPECULAR,   glm::vec3(0.999f),              glm::vec3() };
Material glass =        { Material::DIELECTRIC, glm::vec3(0.999f),              glm::vec3() };
Material light =        { Material::DIFFUSE,    glm::vec3(),                    glm::vec3(0.8f) };

Scene cornellBox =
{
  {
    Sphere { glm::vec3(-0.6f, -0.6f, 1.0f), 0.4f, glass },
    Sphere { glm::vec3(0.0f, -0.6f, 0.0f), 0.4f, mirror },
    Sphere { glm::vec3(0.6f, -0.6f, 1.0f), 0.4f, diffuseGreen },
    Sphere { glm::vec3(0.0f, 1.0 + 9.99f, 0.0f), 10.0f, light },
  },
  {
    Plane { glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, 0.0f, 1.0f), diffuseWhite }, // back
    Plane { glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), diffuseRed }, // left
    Plane { glm::vec3(1.0f,  0.0f, 0.0f), glm::vec3(-1.0f, 0.0f, 0.0f), diffuseBlue }, // right
    Plane { glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), diffuseWhite }, // bottom
    Plane { glm::vec3(0.0f,  1.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f), diffuseWhite }, // top
  }
};

// Intersect a list of objects, returning hit objects closer than outT as hitObj and updating outT
// Make sure to initialise outT first (e.g. to FLT_MAX or the max distance), and the hitObj to null
// so that you can see if there was an intersection
template <typename T>
void intersectObjects(const Ray& ray, const std::vector<T>& objs, const Object*& hitObj, float& outT)
{
  for (auto it = objs.begin(); it != objs.end(); ++it)
  {
    const Object& obj = *it;
    float T;
    if (obj.intersect(ray, T) && T > 0.0f && T < outT)
    {
      hitObj = &obj;
      outT = T;
    }
  }
}

// Tests intersections against a scene description and returns
// the hit object and distance down the ray as hitObj and outT
void intersectScene(const Ray& ray, const Scene& scene, const Object*& hitObj, float& outT)
{
  hitObj = nullptr;
  outT = FLT_MAX;

  intersectObjects(ray, scene.spheres, hitObj, outT);
  intersectObjects(ray, scene.planes, hitObj, outT);
}

// Generates a random float from 0..1
float randf()
{
  static std::default_random_engine e;
  static std::uniform_real_distribution<float> dis(0, 1); // rage 0 - 1
  return dis(e);
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
    u = glm::normalize(glm::vec3(n.x, 0.0f, -n.x));
  else
    u = glm::normalize(glm::vec3(0.0f, -n.z, n.y));
  v = glm::cross(n, u);
}

// Random (uniformly distributed) point on a hemisphere
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

// Transforms a vector into a coordinate space obtained from createCoordinateSpace
glm::vec3 transformIntoSpace(const glm::vec3& sample, const glm::vec3& n, const glm::vec3& u, const glm::vec3& v)
{
  return glm::vec3(
    sample.x * v.x + sample.y * n.x + sample.z * u.x,
    sample.x * v.y + sample.y * n.y + sample.z * u.y,
    sample.x * v.z + sample.y * n.z + sample.z * u.z
  );
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
  if (nearest->mat.type == Material::DIFFUSE)
  {
    // Find coordinate space perpendicular to our normal
    glm::vec3 u, v;
    createCoordinateSpace(hitInfo.nrm, u, v);

    // Create a random sample about a hemisphere
    glm::vec3 sampleDir = uniformSampleHemisphere();

    // Transform it into our normal's coordinate space
    sampleDir = glm::normalize(transformIntoSpace(sampleDir, hitInfo.nrm, u, v));

    Ray sampleRay = Ray { hitInfo.pos, sampleDir };

    glm::vec3 sample = radiance(sampleRay, scene, depth + 1, maxDepth);

    return nearest->mat.emissive + col * sample;
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
      if (rayAngle < criticalAngle)
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
void trace(const Ray& ray, const Scene& scene, int spp, glm::vec4& rgba)
{
  const float avg = (1.0f / spp);

  float pdf = INV_2PI;

  glm::vec3 col;
  for (int i = 0; i < spp; ++i)
    col += avg * radiance(ray, scene, 0, 10) / pdf;

  rgba = glm::vec4(col, 1.0f);
}

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
  glm::vec3 camPos = glm::vec3(0.0f, 0.0f, 5.0f); // The camera position
  float fov = 35.0f; // The fov in degrees
  const int spp = 8; // The number of samples to take per pixel

  // Housekeeping for the progress meter
  float rowPercent = 0.0f;
  float rowPercentStep = 100.0f / height;
  float progressStep = 2;
  float lastProgressVal = 0;

  // Print first part of progress display
  printf("progress: [");

  // Iterate over each pixel and trace scene
  uint8_t* pix = pixels.data();
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      // Construct ray based on fov and pixel
      float px = (2.0f * ((x + 0.5f) / width) - 1.0f) * tan(fov / 2.0f * PI / 180.0f) * aspect;
      float py = (1.0f - 2.0f * ((y + 0.5f) / height)) * tan(fov / 2.0f * PI / 180.0f);
      Ray ray = Ray{ camPos, glm::normalize(glm::vec3(px, py, -1.0f)) };

      // Trace ray
      glm::vec4 rgba;
      trace(ray, sceneToRender, spp, rgba);

      // Clamp result
      rgba = glm::clamp(rgba, glm::vec4(0.0f), glm::vec4(1.0f));

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