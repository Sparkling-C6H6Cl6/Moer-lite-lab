#pragma once

#include "Medium.h"

class BeerslawMedium : public Medium
{
public:
  BeerslawMedium(const Spectrum& density, std::shared_ptr<PhaseFunction> phase);
  BeerslawMedium(const Json& json);

  virtual bool sampleDistance(MediumSampleRecord* rec,
                              const Ray& ray,
                              const Intersection& its,
                              Vector2f sample) const override;

  virtual Spectrum evalTransmittance(Point3f src, Point3f dest) const override;

private:
  Spectrum mDensity;
};

