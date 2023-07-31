#include "HGPhase.h"


PhaseEvalResult
HGPhase::evalPhase(const Vector3f& wo, const Vector3f& wi) const
{
  float phaseValue = phaseHG(dot(normalize(wo), normalize(wi)), g);
  return {
    phaseValue,
    phaseValue
  };
}

PhaseSampleResult
HGPhase::samplePhase(const Vector3f& wo, const Vector2f& sample) const
{
  float sqr = (1.f - g * g) / (1.f + g - 2.f * g * sample[0]);
  float cosTheta = -(1.f + g * g - sqr * sqr) / (2.f * g);
  float sinTheta = std::sqrt(std::max(0.f, 1.f - cosTheta * cosTheta));
  float phi = 2.f * PI * sample[1];

  Vector3f v1 = (std::abs(wo[0]) > std::abs(wo[1]))
    ? normalize(Vector3f{-wo[2], 0, wo[0]})
    : normalize(Vector3f{0, wo[2], -wo[1]});
  Vector3f v2 = cross(wo, v1);
  Vector3f wi = v1 * sinTheta * fm::cos(phi) + v2 * sinTheta * fm::sin(phi) + normalize(wo) * cosTheta;
  float phaseValue = phaseHG(dot(normalize(wo), wi), g);
  return {
    wi,
    phaseValue,
    phaseValue
  };
}

REGISTER_CLASS(HGPhase, "anisotropic")
