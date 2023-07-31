#pragma once

#include "FunctionLayer/Integrator/Integrator.h"
#include "FunctionLayer/Medium/Medium.h"

class VolumePathIntegrator : public Integrator
{
public:
  VolumePathIntegrator() = default;
  VolumePathIntegrator(const Json& json)
    : maxDepth(fetchRequired<uint32_t>(json, "maxDepth")) {}

  virtual ~VolumePathIntegrator() = default;

  virtual Spectrum li(Ray& ray, const Scene& scene,
                      std::shared_ptr<Sampler> sampler) const override;

private:
  struct PathIntegratorRecord
  {
    Vector3f wi;
    Spectrum f;
    double pdf;
    bool isDelta = false;
  };

  static PathIntegratorRecord
    VolumePathIntegrator::evalEmittance(const Scene& scene,
                                        std::optional<Intersection> intersectionOpt,
                                        const Ray& ray);

  static Spectrum evalTransmittance(const Scene& scene,
                                    const Intersection& intersection,
                                    Point3f pointOnLight);

  static Intersection fulfillScatteringIntersection(std::shared_ptr<Medium> medium,
                                                    Point3f position,
                                                    Vector3f normal);

  static PathIntegratorRecord sampleDirectLighting(const Scene& scene,
                                                   std::shared_ptr<Sampler> sampler,
                                                   const Intersection& intersection,
                                                   const Ray& ray);

  static PathIntegratorRecord evalScatter(const Intersection& intersection,
                                          const Ray& ray,
                                          Vector3f wi);

  static PathIntegratorRecord sampleScatter(const Intersection& intersection,
                                            std::shared_ptr<Sampler> sampler,
                                            const Ray& ray);

  static std::shared_ptr<Medium> getTargetMedium(const Intersection& intersection,
                                                 Vector3f wi);

  static std::pair<std::optional<Intersection>, Spectrum>
    intersectIgnoreSurface(const Scene& scene,
                           const Ray& ray,
                           std::shared_ptr<Medium> medium);

  uint32_t maxDepth;
};
