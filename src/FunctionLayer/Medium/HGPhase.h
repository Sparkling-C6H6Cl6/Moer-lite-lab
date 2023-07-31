#pragma once

#include "Phase.h"


class HGPhase : public PhaseFunction
{
  inline float phaseHG(float cosTheta, float g) const {
    float denom = 1.f + g * g + 2.f * g * cosTheta;
    return INV_PI * .25f * (1 - g * g) / (denom * fm::sqrt(denom));
  }
public:
  HGPhase(float _g) : g(_g) {}
  HGPhase(const Json& json) :
    g(std::clamp(fetchOptional(json, "g", 0.f), -1.f, 1.f))
  {}

  virtual PhaseEvalResult evalPhase(const Vector3f& wo, const Vector3f& wi) const override;

  virtual PhaseSampleResult samplePhase(const Vector3f& wo, const Vector2f& sample) const override;

private:
  const float g;
};
