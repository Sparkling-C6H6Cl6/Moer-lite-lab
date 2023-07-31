﻿#pragma once
#include <CoreLayer/Math/Math.h>

inline Vector3f squareToUniformHemisphere(Vector2f sample) {
  float y = 1 - 2 * sample[0];
  float r = fm::sqrt(std::max(.0f, 1.f - y * y));
  float phi = 2 * PI * sample[1];
  Vector3f dir{r * fm::sin(phi), std::abs(y), r * fm::cos(phi)};
  return normalize(dir);
}

inline float squareToUniformHemispherePdf(Vector3f v) {
  return v[1] >= .0f ? INV_PI * .5f : .0f;
}

inline Vector3f squareToCosineHemisphere(Vector2f sample) {
  float phi = 2 * PI * sample[0], theta = fm::acos(std::sqrt(sample[1]));
  return Vector3f{fm::sin(theta) * fm::sin(phi), fm::cos(theta),
                  fm::sin(theta) * fm::cos(phi)};
}

inline float squareToCosineHemispherePdf(Vector3f v) {
  return (v[1] > .0f) ? v[1] * INV_PI : .0f;
}

inline Vector3f squareToUniformSphere(Vector2f sample) {
  float y = 1 - 2 * sample[0];
  float r = fm::sqrt(std::max((float)0, (float)1 - y * y));
  float phi = 2 * PI * sample[1];
  return {r * fm::sin(phi), y, r * fm::cos(phi)};
}

inline float squareToUniformSpherePdf(Vector3f v) {
  return 0.25f * INV_PI;
}
