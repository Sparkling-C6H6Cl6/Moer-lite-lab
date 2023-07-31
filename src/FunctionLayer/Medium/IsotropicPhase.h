#pragma once

#include "Phase.h"
#include "CoreLayer/ColorSpace/Spectrum.h"


class IsotropicPhase : public PhaseFunction
{
public:
  IsotropicPhase() = default;
  IsotropicPhase(const Json& json) {}

  virtual PhaseEvalResult evalPhase(const Vector3f& wo, const Vector3f& wi) const override;

  virtual PhaseSampleResult samplePhase(const Vector3f& wo, const Vector2f& sample) const override;
};
