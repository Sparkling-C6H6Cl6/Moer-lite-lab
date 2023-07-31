#include "VolumePathIntegrator.h"
#include <FunctionLayer/Material/BxDF/Warp.h>
#include <FunctionLayer/Material/Matte.h>


Spectrum
VolumePathIntegrator::li(Ray& ray, const Scene& scene,
                         std::shared_ptr<Sampler> sampler) const
{
  Spectrum spectrum(.0f);
  Spectrum throughput(1.f);
  auto depth = 0u;

  // 照相机可能位于介质中，但是由于底层代码，Shape无法访问Camera的位置，因此暂时不作处理

  auto intersectionOpt = scene.rayIntersect(ray);
  std::shared_ptr<Medium> currentMedium = nullptr;
  Intersection mediumScatteringIntersection;
  MediumSampleRecord record;

  while (true) {

    // eval emittance

    // 光线未击中任何物体
    if (!intersectionOpt.has_value()) {
      for (auto light : scene.infiniteLights) {
        spectrum += throughput * light->evaluateEmission(ray);
      }
      return spectrum;
    }

    Intersection intersection = intersectionOpt.value();
    computeRayDifferentials(&intersection, ray);

    if (depth == 0) {
      // 击中光源
      if (auto light = intersection.shape->light; light) {
        spectrum += light->evaluateEmission(intersection, -ray.direction);
      }
    }

    ++depth;
    if (depth >= maxDepth) {
      break;
    }

    if (currentMedium != nullptr &&
        currentMedium->sampleDistance(&record, ray, intersection, sampler->next2D())) {

      // 光线正位于介质中
      throughput *= record.transmittance * record.sigmaS / record.pdf;
      mediumScatteringIntersection = fulfillScatteringIntersection(currentMedium, record.scatterPoint, ray.direction);

      PathIntegratorRecord lightSampleRecord = sampleDirectLighting(scene, sampler, mediumScatteringIntersection, ray);
      PathIntegratorRecord scatterEvalRecord = evalScatter(mediumScatteringIntersection, ray, lightSampleRecord.wi);
      if (lightSampleRecord.isDelta) {
        spectrum += throughput * lightSampleRecord.f * scatterEvalRecord.f / lightSampleRecord.pdf;
      } else {
        spectrum += throughput * lightSampleRecord.f * scatterEvalRecord.f / (lightSampleRecord.pdf + scatterEvalRecord.pdf);
      }

      // Phase Sample
      PathIntegratorRecord sampleScatterRecord = sampleScatter(mediumScatteringIntersection, sampler, ray);
      if (sampleScatterRecord.f.isZero()) {
        break;
      }
      throughput *= sampleScatterRecord.f / sampleScatterRecord.pdf;

      ray = Ray{mediumScatteringIntersection.position + sampleScatterRecord.wi * 1e-4f,
        sampleScatterRecord.wi};
      intersectionOpt = scene.rayIntersect(ray);

      auto [intersectionOpt, transmittance] = intersectIgnoreSurface(scene, ray, currentMedium);
      PathIntegratorRecord lightEvalRecord = evalEmittance(scene, intersectionOpt, ray);
      if (!lightEvalRecord.f.isZero()) {
        spectrum += throughput * transmittance * lightEvalRecord.f;
      }
    } else {
      if (currentMedium != nullptr) {
        throughput *= record.transmittance / record.pdf;
        currentMedium = nullptr;
      }

      // 光线击中空材质
      if (intersection.shape->material->computeBSDF(intersection) == nullptr) {
        depth--;
        currentMedium = getTargetMedium(intersection, ray.direction);
        ray = Ray{intersection.position + ray.direction * 1e-4f, ray.direction};
        intersectionOpt = scene.rayIntersect(ray);
        continue;
      }

      PathIntegratorRecord lightSampleRecord = sampleDirectLighting(scene, sampler, intersection, ray);
      PathIntegratorRecord scatterEvalRecord = evalScatter(intersection, ray, lightSampleRecord.wi);
      if (lightSampleRecord.isDelta) {
        spectrum += throughput * lightSampleRecord.f * scatterEvalRecord.f / lightSampleRecord.pdf;
      } else {
        spectrum += throughput * lightSampleRecord.f * scatterEvalRecord.f / (lightSampleRecord.pdf + scatterEvalRecord.pdf);
      }

      // Phase Sample
      PathIntegratorRecord sampleScatterRecord = sampleScatter(intersection, sampler, ray);
      if (sampleScatterRecord.f.isZero()) {
        break;
      }
      throughput *= sampleScatterRecord.f / sampleScatterRecord.pdf;
      currentMedium = getTargetMedium(intersection, sampleScatterRecord.wi);

      ray = Ray{intersection.position + sampleScatterRecord.wi * 1e-4f,
        sampleScatterRecord.wi};
      intersectionOpt = scene.rayIntersect(ray);

      auto [intersectionOpt, transmittance] = intersectIgnoreSurface(scene, ray, currentMedium);
      PathIntegratorRecord lightEvalRecord = evalEmittance(scene, intersectionOpt, ray);
      if (!lightEvalRecord.f.isZero()) {
        spectrum += throughput * transmittance * lightEvalRecord.f;
      }
    }
  }
  return spectrum;
}


VolumePathIntegrator::PathIntegratorRecord
VolumePathIntegrator::evalEmittance(const Scene& scene,
                                    std::optional<Intersection> intersectionOpt,
                                    const Ray& ray)
{
  Vector3f wo = -ray.direction;
  Spectrum spectrum(0.f);
  if (!intersectionOpt.has_value())
  {
    float pdf = 0.f;
    for (auto light : scene.infiniteLights) {
      spectrum += light->evaluateEmission(ray);

    }
  } else if (auto light = intersectionOpt.value().shape->light; light) {
    Intersection intersection = intersectionOpt.value();
    Vector3f n = intersection.normal;
    spectrum = light->evaluateEmission(intersection, ray.direction);
  }
  return {
    ray.direction,
    spectrum,
    0.f,
    false
  };
}


Spectrum
VolumePathIntegrator::evalTransmittance(const Scene& scene,
                                        const Intersection& intersection,
                                        Point3f pointOnLight)
{
  Vector3f direction = pointOnLight - intersection.position;
  float maxDistance = direction.length();
  direction = normalize(direction); // -wi

  Spectrum transmittance(1.f);
  Ray ray{intersection.position + direction * 1e-4f, direction};
  auto medium = getTargetMedium(intersection, -direction);
  auto testRayIntersectionOpt = scene.rayIntersect(ray);

  Point3f target = pointOnLight;
  Point3f lastScatteringPoint = intersection.position;

  // 光源沿直线无限传播，逐个计算透射率
  while (true) {

    // 特殊情况：光线传播无限远或无限光源
    if (!testRayIntersectionOpt.has_value()) {
      if (medium != nullptr) {
        transmittance *= medium->evalTransmittance(target, lastScatteringPoint);
      }
      return transmittance;
    }

    const Intersection& testRayIntersection = testRayIntersectionOpt.value();

    // 特殊情况：光线击中光源
    if ((testRayIntersection.position - intersection.position).length()
        >= maxDistance - 1e-4f) {
      if (medium != nullptr) {
        transmittance *= medium->evalTransmittance(target, lastScatteringPoint);
      }
      return transmittance;
    }

    // 特殊情况：光线击中非空材质
    if (testRayIntersection.shape->material != nullptr &&
        testRayIntersection.shape->material->computeBSDF(testRayIntersection) != nullptr) {
      return 0.f;
    }

    if (medium != nullptr) {
      transmittance *= medium->evalTransmittance(testRayIntersection.position, lastScatteringPoint);
    }

    medium = getTargetMedium(testRayIntersection, direction);
    ray.origin = testRayIntersection.position + direction * 1e-4f;
    lastScatteringPoint = testRayIntersection.position;
    testRayIntersectionOpt = scene.rayIntersect(ray);
  }
}


Intersection
VolumePathIntegrator::fulfillScatteringIntersection(std::shared_ptr<Medium> medium,
                                                    Point3f position,
                                                    Vector3f normal)
{
  Intersection scatteringPoint;
  scatteringPoint.position = position;
  scatteringPoint.normal = -normal;
  scatteringPoint.scatteringMedium = medium;
  return scatteringPoint;
}


VolumePathIntegrator::PathIntegratorRecord
VolumePathIntegrator::sampleDirectLighting(const Scene& scene,
                                           std::shared_ptr<Sampler> sampler,
                                           const Intersection& intersection,
                                           const Ray& ray)
{
  float pdfLight = .0f;
  auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
  LightSampleResult result = light->sample(intersection, sampler->next2D());
  float pdfDirect = result.pdf * pdfLight;
  Vector3f wi = result.direction;
  Point3f pointOnLight = intersection.position + result.direction * result.distance;
  auto transmittance = evalTransmittance(scene, intersection, pointOnLight);
  return {
    wi,
    transmittance * result.energy,
    pdfDirect,
    result.isDelta
  };
}


VolumePathIntegrator::PathIntegratorRecord
VolumePathIntegrator::evalScatter(const Intersection& intersection,
                                  const Ray& ray,
                                  Vector3f wi)
{
  if (intersection.scatteringMedium == nullptr) {
    auto bsdf = intersection.shape->material->computeBSDF(intersection);
    Vector3f n = intersection.normal;
    double wiDotN = fm::abs(dot(n, wi));
    Vector3f wo = -ray.direction;
    return {
      wi,
      bsdf->f(wo, wi) * wiDotN,
      bsdf->pdf(wo, wi),
      false
    };
  } else {
    auto medium = intersection.scatteringMedium;
    Vector3f wo = -ray.direction;
    PhaseEvalResult result = medium->evalPhase(wo, wi);
    return {
      wi,
      result.value,   // no cosine term
      result.pdf,
      false
    };
  }
}


VolumePathIntegrator::PathIntegratorRecord
VolumePathIntegrator::sampleScatter(const Intersection& intersection,
                                    std::shared_ptr<Sampler> sampler,
                                    const Ray& ray)
{
  if (intersection.scatteringMedium == nullptr) {
    Vector3f wo = -ray.direction;
    auto bsdf = intersection.shape->material->computeBSDF(intersection);
    Vector3f n = intersection.normal;
    BSDFSampleResult result = bsdf->sample(wo, sampler->next2D());
    float pdf = result.pdf;
    Vector3f wi = result.wi;
    double wiDotN = fm::abs(dot(wi, n));
    return {
      wi,
      result.weight * wiDotN,
      pdf,
      result.type == BSDFType::Specular
    };
  } else {
    Vector3f wo = -ray.direction;
    auto medium = intersection.scatteringMedium;
    Point3f scatteringPoint = intersection.position;
    PhaseSampleResult result = medium->samplePhase(wo, sampler->next2D());
    Spectrum phaseValue = result.value;
    double pdf = result.pdf;
    Vector3f wi = result.wi;
    return {
      wi,
      phaseValue,
      pdf,
      false
    };
  }
}


std::shared_ptr<Medium>
VolumePathIntegrator::getTargetMedium(const Intersection& intersection,
                                      Vector3f wi)
{
  bool scatterToOutSide = dot(intersection.normal, wi) > 0;
  if (intersection.scatteringMedium == nullptr) {
    if (scatterToOutSide) {
      return intersection.shape->material->getOutsideMedium();
    } else {
      return intersection.shape->material->getInsideMedium();
    }
  } else {
    return nullptr;
  }
}


std::pair<std::optional<Intersection>, Spectrum>
VolumePathIntegrator::intersectIgnoreSurface(const Scene& scene,
                                             const Ray& ray,
                                             std::shared_ptr<Medium> medium)
{
  Vector3f direction = ray.direction;
  Spectrum transmittance(0.f);
  Ray marchRay{ray.origin + direction * 1e-4f, direction};
  auto currentMedium = medium;
  Point3f lastScatteringPoint = ray.origin;
  auto testRayIntersectionOpt = scene.rayIntersect(marchRay);

  // calculate the transmittance of last segment from lastScatteringPoint to testRayItsOpt.
  while (true) {

    // corner case: infinite medium or infinite light source.
    if (!testRayIntersectionOpt.has_value()) {
      if (currentMedium != nullptr) {
        transmittance = Spectrum(0.0);
      }
      return {testRayIntersectionOpt, transmittance};
    }

    Intersection testRayIntersection = testRayIntersectionOpt.value();

    // corner case: non-null surface
    if (testRayIntersection.scatteringMedium == nullptr) {
      if (testRayIntersection.shape->material->computeBSDF(testRayIntersection) == nullptr) {
        if (currentMedium != nullptr)
          transmittance *= currentMedium->evalTransmittance(testRayIntersection.position, lastScatteringPoint);
        return {testRayIntersectionOpt, transmittance};
      }
    }

    // hit a null surface, calculate tr
    if (currentMedium != nullptr) {
      transmittance *= currentMedium->evalTransmittance(testRayIntersection.position, lastScatteringPoint);
    }

    // update medium
    currentMedium = getTargetMedium(testRayIntersection, direction);

    // update ray and intersection point.
    marchRay.origin = testRayIntersection.position + direction * 1e-4f;
    lastScatteringPoint = testRayIntersection.position;
    testRayIntersectionOpt = scene.rayIntersect(marchRay);
  }
  return {testRayIntersectionOpt, transmittance};
}


REGISTER_CLASS(VolumePathIntegrator, "volumePath")
