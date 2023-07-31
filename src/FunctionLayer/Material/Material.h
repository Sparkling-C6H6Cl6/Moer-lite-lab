#pragma once
#include "./BxDF/BSDF.h"
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Shape/Intersection.h>
#include <FunctionLayer/Texture/NormalTexture.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/JsonUtil.h>

#include "FunctionLayer/Medium/Medium.h"
#include "FunctionLayer/Medium/Beerslaw.h"
#include "FunctionLayer/Medium/IsotropicPhase.h"


class Material {
public:
  Material() {
    // donothing
  }

  Material(const Json& json) {
    if (json.contains("normalmap")) {
      normalMap = std::make_shared<NormalTexture>(json["normalmap"]);
    }
    twoSideShading = fetchOptional(json, "twoSideShading", true);
    if (json.contains("insideMedium")) {
      insideMedium = Factory::construct_class<Medium>(json["insideMedium"]);
    } else {
      insideMedium = nullptr;
    }
    if (json.contains("outsideMedium")) {
      outsideMedium = Factory::construct_class<Medium>(json["outsideMedium"]);
    } else {
      outsideMedium = nullptr;
    }
  }

  virtual std::shared_ptr<BSDF> computeBSDF(const Intersection& intersection) const = 0;

  void computeShadingGeometry(const Intersection& intersection,
                              Vector3f* normal, Vector3f* tangent,
                              Vector3f* bitangent) const;

  inline std::shared_ptr<Medium> getInsideMedium() const {
    return insideMedium;
  }

  inline std::shared_ptr<Medium> getOutsideMedium() const {
    return outsideMedium;
  }

  inline void setInsideMedium(std::shared_ptr<Medium> _insideMedium) {
    insideMedium = _insideMedium;
  }

  inline void setOutsideMedium(std::shared_ptr<Medium> _outsideMedium) {
    outsideMedium = _outsideMedium;
  }

protected:
  std::shared_ptr<NormalTexture> normalMap;
  bool twoSideShading;
  std::shared_ptr<Medium> insideMedium;
  std::shared_ptr<Medium> outsideMedium;
};