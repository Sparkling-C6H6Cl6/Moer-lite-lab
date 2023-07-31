#pragma once

#include "Material.h"

class NullMaterial : Material
{
public:
  NullMaterial();
  NullMaterial(const Json& json);
  virtual std::shared_ptr<BSDF> computeBSDF(const Intersection& intersection) const override;
};
