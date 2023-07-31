#pragma once

#include "CoreLayer/Math/Geometry.h"
#include "CoreLayer/Math/Math.h"
#include "ResourceLayer/Factory.h"

struct PhaseEvalResult
{
  float value;
  float pdf;
};

struct PhaseSampleResult
{
  Vector3f wi;
  float value;
  float pdf;
};

class PhaseFunction
{
public:

  /**
   * @brief Given the direction of wi/wo (both in world),
   *        evaluate the phase function
   *
   * @param wo The direction towards light/next intersection
   * @param wi The direction towards camera/previous intersection
   */
  virtual PhaseEvalResult evalPhase(const Vector3f& wo, const Vector3f& wi) const = 0;

  /**
   * @brief Given the direction of wo (in world),
   *        sample a scatter direction
   *
   * @param wo
   * @param scatterPoint
   */
  virtual PhaseSampleResult samplePhase(const Vector3f& wo, const Vector2f& sample) const = 0;
};
