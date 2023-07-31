#include "NullMaterial.h"

NullMaterial::NullMaterial() {}

NullMaterial::NullMaterial(const Json& json) : Material(json) {}

std::shared_ptr<BSDF>
NullMaterial::computeBSDF(const Intersection& intersection) const { return nullptr; }

REGISTER_CLASS(NullMaterial, "null")
