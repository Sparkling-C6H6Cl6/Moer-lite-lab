#pragma once
#include "BSDF.h"
#include "Warp.h"

class OrenNayarBSDF : public BSDF {

  template <
    typename T,
    typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
  constexpr static T _sq(T x) {
    return x * x;
  }

public:
  OrenNayarBSDF(const Vector3f& _normal, const Vector3f& _tangent,
                const Vector3f& _bitangent, Spectrum _albedo, float _sigma)
    : BSDF(_normal, _tangent, _bitangent), albedo(_albedo), sigma(_sigma) {}

  virtual Spectrum f(const Vector3f& wo, const Vector3f& wi) const override {

    // TODO
    // 1. 转换坐标系到局部坐标

    Vector3f woLocal = toLocal(wo);
    Vector3f wiLocal = toLocal(wi);

    // 2. 计算 A, B, \alpha, \beta（可以直接求\sin\alpha,\tan\beta）,
    // \cos(\phi_i-\phi_o)

    float paramA = 1.f - (.5f * _sq(sigma)) / (_sq(sigma) + .33f);
    float paramB = (.45f * _sq(sigma)) / (_sq(sigma) + .09f);

    float cosThetaO = clamp(dot(wo, Vector3f{0.f, 1.f, 0.f}), 0.f, 1.f);
    float cosThetaI = clamp(dot(wi, Vector3f{0.f, 1.f, 0.f}), 0.f, 1.f);

    float cosAlpha = std::min(cosThetaI, cosThetaO);
    float cosBeta = std::max(cosThetaI, cosThetaO);
    float sinAlpha = fm::sqrt(1.f - _sq(cosAlpha));
    float tanBeta = std::numeric_limits<float>::max();
    if (cosBeta != 0.f) {
      tanBeta = fm::sqrt(1.f - _sq(cosBeta)) / cosBeta;
    }

    float cosGamma = dot(normalize(wo * Vector3f{1.f, 0.f, 1.f}),
                         normalize(wi * Vector3f{1.f, 0.f, 1.f}));

    // 3. return Oren-Nayar brdf

    return albedo / PI *
      (paramA + paramB * std::max(0.f, cosGamma) * sinAlpha * tanBeta);

  }

  virtual BSDFSampleResult sample(const Vector3f& wo,
                                  const Vector2f& sample) const override {
    Vector3f wi = squareToCosineHemisphere(sample);
    float pdf = squareToCosineHemispherePdf(wi);
    return {albedo, toWorld(wi), pdf, BSDFType::Diffuse};
  }

private:
  Spectrum albedo;
  float sigma;
};