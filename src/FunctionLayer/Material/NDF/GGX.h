#pragma once

#include "NDF.h"

class GGXDistribution : public NDF {

  template <
    typename T,
    typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr static T _sq(T x) {
    return x * x;
  }

  float getG1(const Vector3f& w, const Vector2f& alpha) const {

    float cosTheta = w[1];
    float tanTheta = fm::sqrt(1.f - _sq(cosTheta)) / cosTheta;

    Vector2f vecPhi = Vector2f{w[0], w[2]};
    vecPhi /= vecPhi.len();

    float alphaNum = fm::sqrt(alpha[0] * _sq(vecPhi[0]) + alpha[1] * _sq(vecPhi[1]));

    return 2.f / (1.f + fm::sqrt(1.f + _sq(alphaNum) * _sq(tanTheta)));
  }

public:

  GGXDistribution() noexcept = default;

  virtual ~GGXDistribution() noexcept = default;

  virtual float getD(const Vector3f& whLocal,
                     const Vector2f& alpha) const noexcept override {

    float cosTheta = whLocal[1];
    float tanTheta = fm::sqrt(1.f - _sq(cosTheta)) / cosTheta;

    Vector2f vecPhi = Vector2f(whLocal[0], whLocal[2]);
    vecPhi /= vecPhi.len();

    float alphaNum = fm::sqrt(alpha[0] * _sq(vecPhi[0]) + alpha[1] * _sq(vecPhi[1]));

    return _sq(alphaNum) /
      (PI * _sq(_sq(cosTheta)) * _sq(_sq(alphaNum) + _sq(tanTheta)));
  }

  virtual float getG(const Vector3f& woLocal, const Vector3f& wiLocal,
                     const Vector2f& alpha) const noexcept override {
    return getG1(woLocal, alpha) * getG1(wiLocal, alpha);
  }

  virtual float pdf(const Vector3f& woLocal, const Vector3f& whLocal,
                    const Vector2f& alpha) const noexcept override {
    return getD(whLocal, alpha) * whLocal[1];
  }

  virtual Vector3f sampleWh(const Vector3f& woLocal, const Vector2f& alpha,
                            const Vector2f& sample) const noexcept override {
    float a = alpha[0];
    float tan_theta_2 = a * a * sample[0] / (1.f - sample[0]);
    float phi = sample[1] * 2 * PI;

    float cos_theta = std::sqrt(1.f / (1.f + tan_theta_2));
    float sin_theta = std::sqrt(1.f - cos_theta * cos_theta);
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
  }

};