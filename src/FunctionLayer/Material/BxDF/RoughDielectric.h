#pragma once
#include "../NDF/NDF.h"
#include "BSDF.h"
#include "Warp.h"
#include <random>

class RoughDielectricBSDF : public BSDF {

  template <
    typename T,
    typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr static T _sq(T x) {
    return x * x;
  }

public:
  RoughDielectricBSDF(const Vector3f& _normal, const Vector3f& _tangent,
                      const Vector3f& _bitangent, Spectrum _albedo,
                      Vector2f _alpha, float _eta, NDF* _ndf)
    : BSDF(_normal, _tangent, _bitangent), albedo(_albedo), alpha(_alpha),
    eta(_eta), ndf(_ndf) {}

  virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const override {

    // TODO
    // 1. 转换坐标系到局部坐标

    Vector3f woLocal = normalize(toLocal(wo));
    Vector3f wiLocal = normalize(toLocal(wi));
    Vector3f whLocal = normalize(woLocal + wiLocal);

    // 2. 根据公式计算 Fr, D, G

    float cosThetaO = woLocal[1];
    float cosThetaI = wiLocal[1];
    float etaO = cosThetaI >= 0.f ? eta : (1.f / eta);
    float paramFr = getFr(etaO, dot(woLocal, whLocal));
    float paramD = ndf->getD(whLocal, alpha);
    float paramG = ndf->getG(woLocal, wiLocal, alpha);

    // 3. return albedo * D * G * Fr / (4 * \cos\theta_o);
    // tips:
    // 不考虑多重介质，如果光线从真空射入介质，其eta即配置中填写的eta；
    // 如果光线从介质射出，则eta = 1/eta

    return albedo * paramD * paramG * paramFr / (4.f * cosThetaO);
  }

  virtual float pdf(const Vector3f& wo, const Vector3f& wi) const override {
    Vector3f woLocal = toLocal(wo);
    Vector3f wiLocal = toLocal(wi);
    Vector3f whLocal = normalize(woLocal + wiLocal);
    float woDotWh = dot(woLocal, whLocal);
    if (woDotWh < 0.f || wiLocal[1] < 0.f) {
      float cosThetaT = fm::sqrt(std::max(1.f - _sq(eta) * (1.f - _sq(woDotWh)), 0.f));
      Vector3f glossLocal = (eta * woDotWh - (woDotWh > 0.f ? 1.f : -1.f) * cosThetaT) * whLocal - eta * woLocal;
      return ndf->pdf(woLocal, whLocal, alpha) * (1 - getFr(eta, woDotWh)) *
        std::abs(dot(glossLocal, whLocal) / _sq(dot(woLocal, whLocal) * eta + dot(glossLocal, whLocal)));
    } else {
      return ndf->pdf(woLocal, whLocal, alpha) * 0.25f / woDotWh;
    }
  }

  virtual BSDFSampleResult sample(const Vector3f& wo,
                                  const Vector2f& sample) const override {
    Vector3f woLocal = toLocal(wo);
    Vector3f whLocal = ndf->sampleWh(woLocal, alpha, sample);
    Vector3f wiLocal = normalize(-woLocal + 2.f * dot(woLocal, whLocal) * whLocal);
    float woDotWh = dot(woLocal, whLocal);
    float paramFr = getFr(eta, woDotWh);
    float paramD = ndf->getD(whLocal, alpha);
    bool reflect = std::uniform_real_distribution(0.f, 1.f)(std::default_random_engine()) < paramFr;
    if (reflect) {
      float paramG = ndf->getG(woLocal, wiLocal, alpha);
      return {
        paramD * paramG * paramFr / (4.f * woLocal[1] * wiLocal[1]),
        toWorld(wiLocal),
        ndf->pdf(woLocal, whLocal, alpha) * 0.25f / woDotWh,
        BSDFType::GlossyReflection
      };
    } else {
      float cosThetaT = fm::sqrt(std::max(1.f - _sq(eta) * (1.f - _sq(woDotWh)), 0.f));
      Vector3f glossLocal = (eta * woDotWh - (woDotWh > 0.f ? 1.f : -1.f) * cosThetaT) * whLocal - eta * woLocal;
      float paramG = ndf->getG(woLocal, glossLocal, alpha);
      return {
        paramD * paramG * paramFr / (4.f * woLocal[1] * glossLocal[1]),
        toWorld(glossLocal),
        ndf->pdf(woLocal, whLocal, alpha) * (1 - paramFr) *
          std::abs(dot(glossLocal, whLocal) / _sq(dot(woLocal, whLocal) * eta + dot(glossLocal, whLocal))),
        BSDFType::GlossyTransmission
      };
    }
  }

  float getR0(float etaO) const noexcept {
    return ((etaO - 1.f) * (etaO - 1.f)) / ((etaO + 1.f) * (etaO + 1.f));
  }

  float getFr(float etaO, float cosTheta) const noexcept {
    float r0 = getR0(etaO);
    return r0 + (1.f - r0) * std::pow(1.f - cosTheta, 5.f);
  }

private:
  Spectrum albedo;
  Vector2f alpha;
  float eta;
  NDF* ndf;
};