#pragma once

#include "NDF.h"

class BeckmannDistribution : public NDF {

  template <
    typename T,
    typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr static T _sq(T x) {
    return x * x;
  }

  float getG1(const Vector3f& w, const Vector2f& alpha) const {

    float cosTheta = dot(w, Vector3f{0.f, 1.f, 0.f});
    float tanTheta = fm::sqrt(1.f - _sq(cosTheta)) / cosTheta;

    Vector2f phiVec = Vector2f(w[0], w[2]);
    phiVec /= phiVec.len();

    float alphaNum = phiVec.dot(alpha);
    float a = 1.f / (alphaNum * tanTheta);

    if (a < 1.6f) {
      return (2.181f * _sq(a) + 3.535f * a) /
        (2.577f * _sq(a) + 2.276f * a + 1.f);
    } else {
      return 1.f;
    }
  }

public:

  BeckmannDistribution() noexcept = default;

  virtual ~BeckmannDistribution() noexcept = default;

  virtual float getD(const Vector3f& whLocal,
                     const Vector2f& alpha) const noexcept override {

    float cosTheta = dot(whLocal, Vector3f{0.f, 1.f, 0.f});
    float tanTheta = fm::sqrt(1.f - _sq(cosTheta)) / cosTheta;

    Vector2f phiVec = Vector2f(whLocal[0], whLocal[2]);
    phiVec /= phiVec.len();

    float alphaNum = phiVec.dot(alpha);

    return fm::exp(-_sq(tanTheta) / _sq(alphaNum))
      / (PI * _sq(alphaNum) * _sq(_sq(cosTheta)));
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
    float tan_theta_2 = -std::log(1 - sample[0]) * a * a;
    float phi = sample[1] * 2 * PI;

    float cos_theta = std::sqrt(1.f / (1.f + tan_theta_2));
    float sin_theta = std::sqrt(1.f - cos_theta * cos_theta);
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta};
  }

};