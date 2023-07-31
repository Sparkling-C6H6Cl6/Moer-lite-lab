#pragma once

#include "Medium.h"

class HomogeneousMedium : public Medium
{
public:
  HomogeneousMedium(Spectrum sigmaT, Spectrum albedo, std::shared_ptr<PhaseFunction> phase);
  HomogeneousMedium(const Json& json);

  virtual bool sampleDistance(MediumSampleRecord* rec,
                              const Ray& ray,
                              const Intersection& its,
                              Vector2f sample) const override;

  virtual Spectrum evalTransmittance(Point3f src, Point3f dest) const override;

private:
  Spectrum mSigmaA; // 吸收截面（absorption cross section），光在介质中移动的单位距离中被吸收的概率密度
  Spectrum mSigmaS; // 散射系数（scattering coefficient），光在介质中移动的单位距离中被外散射的概率密度
  Spectrum mSigmaT; // 衰减系数（attenuation coefficient），sigma_t = sigma_a + sigma_s
  Spectrum mAlbedo; // 反照率（albedo），albedo = sigma_s / sigma_t
};
