#include "IsotropicPhase.h"
#include "FunctionLayer/Material/BxDF/Warp.h"


PhaseEvalResult
IsotropicPhase::evalPhase(const Vector3f& wo, const Vector3f& wi) const
{
  return {
    .25f * INV_PI,
    .25f * INV_PI
  };
}

PhaseSampleResult
IsotropicPhase::samplePhase(const Vector3f& wo, const Vector2f& sample) const
{
  return {
    squareToUniformSphere(sample),
    .25f * INV_PI,
    .25f * INV_PI
  };
}

REGISTER_CLASS(IsotropicPhase, "isotropic")
