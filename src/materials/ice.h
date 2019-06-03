#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_ICE_H
#define PBRT_MATERIALS_ICE_H

#include "pbrt.h"
#include "material.h"
#include "sampling.h"
#include "imageio.h"
#include "spectrum.h"

namespace pbrt {

class IceMaterial : public Material {
  public:
    IceMaterial(std::string& mapName) {
      Point2i res;
      std::unique_ptr<RGBSpectrum[]> map = ReadImage(mapName, &res);
      Float* data = new Float[res.x * res.y];
      for (int i = 0; i < res.x * res.y; i++) {
        RGBSpectrum rgb = map[i];
        data[i] = rgb[0];
      }
      distr = new Distribution2D(data, res.x, res.y);
    }
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;
  private:
    Distribution2D* distr;
};

IceMaterial *CreateIceMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_ICE_H
