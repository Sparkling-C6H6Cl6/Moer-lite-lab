#include "Homogeneous.h"


HomogeneousMedium::HomogeneousMedium(Spectrum sigmaT,
                                     Spectrum albedo,
                                     std::shared_ptr<PhaseFunction> phase) :
  Medium(phase),
  mAlbedo(albedo),
  mSigmaT(sigmaT),
  mSigmaS(mSigmaT* mAlbedo),
  mSigmaA(mSigmaT - mSigmaS) {}

HomogeneousMedium::HomogeneousMedium(const Json& json) :
  Medium(json)
{
  if (json.contains("albedo")) {
    if (json["albedo"].is_number_float()) {
      mAlbedo = Spectrum(std::clamp(fetchRequired<float>(json, "albedo"), 0.f, 1.f));
    } else {
      mAlbedo = fetchRequired<Spectrum>(json, "albedo");
      mAlbedo[0] = std::clamp(mAlbedo[0], 0.f, 1.f);
      mAlbedo[1] = std::clamp(mAlbedo[1], 0.f, 1.f);
      mAlbedo[2] = std::clamp(mAlbedo[2], 0.f, 1.f);
    }
  } else {
    mAlbedo = Spectrum(0.1f);
  }
  if (json.contains("sigmaT")) {
    if (json["sigmaT"].is_number_float()) {
      mSigmaT = Spectrum(std::max(fetchRequired<float>(json, "sigmaT"), 0.f));
    } else {
      mSigmaT = fetchRequired<Spectrum>(json, "sigmaT");
      mSigmaT[0] = std::max(mSigmaT[0], 0.f);
      mSigmaT[1] = std::max(mSigmaT[1], 0.f);
      mSigmaT[2] = std::max(mSigmaT[2], 0.f);
    }
  } else {
    mSigmaT = Spectrum(0.8f);
  }
  mSigmaS = mSigmaT * mAlbedo;
  mSigmaA = mSigmaT - mSigmaS;
}

bool
HomogeneousMedium::sampleDistance(MediumSampleRecord* rec,
                                  const Ray& ray,
                                  const Intersection& its,
                                  Vector2f sample) const
{
  // randomly pick a channel/frequency and sample a distance.
  int channelIndex = int(sample[0] * 3.f);
  float dist = -fm::log(1.f - sample[1]) / mSigmaT[channelIndex];

  if (dist < its.distance) {
    // sampled a scattering point inside the medium.
    rec->marchLength = dist;
    rec->scatterPoint = ray.at(dist);
    rec->sigmaA = mSigmaA;
    rec->sigmaS = mSigmaS;
    rec->transmittance = evalTransmittance(ray.origin, rec->scatterPoint);
    // calculate pdf
    rec->pdf = (mSigmaT[0] * fm::exp(-mSigmaT[0] * dist)
                + mSigmaT[1] * fm::exp(-mSigmaT[1] * dist)
                + mSigmaT[2] * fm::exp(-mSigmaT[2] * dist)) / 3.f;
    return true;
  } else {
    // sampled a point on object boundary (surface).
    rec->marchLength = its.distance;
    rec->transmittance = evalTransmittance(ray.origin, its.position);
    // calculate discrete probility (instead of continuous probability density)
    rec->pdf = (fm::exp(-mSigmaT[0] * its.distance)
                + fm::exp(-mSigmaT[1] * its.distance)
                + fm::exp(-mSigmaT[2] * its.distance)) / 3.f;
    return false;
  }
  // Incomprehensible strategy that discrete probabilities and continuous probabilities share the same responsibility.
}

Spectrum
HomogeneousMedium::evalTransmittance(Point3f src, Point3f dest) const
{
  float distance = (dest - src).length();
  return Spectrum(fm::exp(-mSigmaT[0] * distance),
                  fm::exp(-mSigmaT[1] * distance),
                  fm::exp(-mSigmaT[2] * distance));
}

REGISTER_CLASS(HomogeneousMedium, "homogeneous")
