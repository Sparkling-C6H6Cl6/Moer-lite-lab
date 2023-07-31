#include "Beerslaw.h"


BeerslawMedium::BeerslawMedium(const Spectrum& density, std::shared_ptr<PhaseFunction> phase) :
  Medium(phase),
  mDensity(density) {}

BeerslawMedium::BeerslawMedium(const Json& json) :
  Medium(json)
{
  if (json.contains("intensity")) {
    if (json["intensity"].is_number_float()) {
      mDensity = Spectrum(std::max(fetchRequired<float>(json, "intensity"), 0.f));
    } else {
      mDensity = fetchRequired<Spectrum>(json, "intensity");
      mDensity[0] = std::max(mDensity[0], 0.f);
      mDensity[1] = std::max(mDensity[1], 0.f);
      mDensity[2] = std::max(mDensity[2], 0.f);
    }
  } else {
    mDensity = Spectrum(1.f);
  }
}

bool
BeerslawMedium::sampleDistance(MediumSampleRecord* rec,
                               const Ray& ray,
                               const Intersection& its,
                               Vector2f sample) const
{
  // No scattering, only absorbtion
  rec->marchLength = its.distance;
  rec->pdf = 1.f;
  rec->transmittance = evalTransmittance(ray.origin, its.position);
  return false;
}

Spectrum
BeerslawMedium::evalTransmittance(Point3f src, Point3f dest) const
{
  float distance = (dest - src).length();
  return Spectrum(fm::exp(mDensity[0] * -distance),
                  fm::exp(mDensity[1] * -distance),
                  fm::exp(mDensity[2] * -distance));
}

REGISTER_CLASS(BeerslawMedium, "beerslaw")
