#pragma once

#include "Phase.h"
#include "CoreLayer/ColorSpace/Spectrum.h"
#include "FunctionLayer/Ray/Ray.h"
#include "FunctionLayer/Shape/Intersection.h"

struct MediumSampleRecord {
  float marchLength;
  float pdf;
  Point3f scatterPoint;
  Vector3f wi;
  Spectrum transmittance;
  Spectrum sigmaA;
  Spectrum sigmaS;
};

class Medium
{
public:
  Medium() = default;
  Medium(const Json& json) {
    mPhase = Factory::construct_class<PhaseFunction>(json["phase"]);
  }
  Medium(std::shared_ptr<PhaseFunction> phase) : mPhase(phase) {}

  /**
   * @brief Sample a distance, the ray will transport in media without collision until reach the distance
   * @return
   * - true, if the ray will endure a collision in medium
   * - false, if the ray will pass through the media without collision
   */
  virtual bool sampleDistance(MediumSampleRecord* rec,
                              const Ray& ray,
                              const Intersection& its,
                              Vector2f sample) const = 0;

  virtual Spectrum evalTransmittance(Point3f src, Point3f dest) const = 0;

  PhaseEvalResult evalPhase(Vector3f wo, Vector3f wi) const {
    return mPhase->evalPhase(wo, wi);
  }

  PhaseSampleResult samplePhase(Vector3f wo, Vector2f sample) const {
    return mPhase->samplePhase(wo, sample);
  }

protected:
  std::shared_ptr<PhaseFunction> mPhase;
};
