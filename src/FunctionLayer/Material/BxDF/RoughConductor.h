#pragma once
#include "../NDF/NDF.h"
#include "BSDF.h"
#include "Warp.h"

class RoughConductorBSDF : public BSDF {

  template <
    typename T,
    typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr static T _sq(T x) {
    return x * x;
  }

public:
  RoughConductorBSDF(const Vector3f& _normal, const Vector3f& _tangent,
                     const Vector3f& _bitangent, Spectrum _albedo,
                     Vector2f _alpha, Vector3f _eta, Vector3f _k, NDF* _ndf)
    : BSDF(_normal, _tangent, _bitangent), albedo(_albedo), alpha(_alpha),
    eta(_eta), k(_k), ndf(_ndf) {}

  virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const override {

    // TODO
    // 1. 转换坐标系到局部坐标

    Vector3f woLocal = normalize(toLocal(wo));
    Vector3f wiLocal = normalize(toLocal(wi));
    Vector3f whLocal = normalize(woLocal + wiLocal);

    // 2. 根据公式计算 Fr, D, G

    float cosThetaO = woLocal[1];
    float cosThetaI = wiLocal[1];

    if (cosThetaO == 0.f || cosThetaI == 0.f || whLocal.isZero()) {
      return 0.f;
    }

    Vector3f paramFr = getFr(dot(woLocal, whLocal));
    float paramD = ndf->getD(whLocal, alpha);
    float paramG = ndf->getG(woLocal, wiLocal, alpha);

    // 3. return albedo * D * G * Fr / (4 * \cos\theta_o);
    // tips: brdf
    // 中分母的\cos\theta_i项被渲染方程中的cos项消去，因此分母只有4*\cos\theta_o

    return albedo * paramD * paramG * paramFr / (4.f * cosThetaO);
  }

  virtual BSDFSampleResult sample(const Vector3f& wo,
                                  const Vector2f& sample) const override {
    Vector3f woLocal = toLocal(wo);
    Vector3f whLocal = ndf->sampleWh(woLocal, alpha, sample);
    Vector3f wiLocal = normalize(-woLocal + 2.f * dot(woLocal, whLocal) * whLocal);
    float woDotWh = dot(woLocal, whLocal);
    if (woDotWh < 0.f || wiLocal[1] < 0.f) {
      return {};
    }
    float pdf = ndf->pdf(woLocal, whLocal, alpha) * 0.25f / woDotWh;
    Vector3f paramFr = getFr(woDotWh);
    float paramD = ndf->getD(whLocal, alpha);
    float paramG = ndf->getG(woLocal, wiLocal, alpha);
    return {
      paramD * paramG * paramFr / (4.f * woLocal[1] * wiLocal[1]),
      toWorld(wiLocal),
      pdf,
      BSDFType::GlossyReflection
    };
  }

  Vector3f getR0() const noexcept {
    return ((eta - 1.f) * (eta - 1.f) + k * k) /
      ((eta + 1.f) * (eta + 1.f) + k * k);
  }

  Vector3f getFr(float cosTheta) const noexcept {
    Vector3f r0 = getR0();
    return r0 + (Vector3f(1.f) - r0) * std::pow(1.f - cosTheta, 5.f);
  }

private:
  Spectrum albedo;
  Vector2f alpha;
  Vector3f eta;
  Vector3f k;
  NDF* ndf;
};