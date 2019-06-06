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
      Float* radi = new Float[res.x * res.y];

      for (int i = 0; i < res.x * res.y; i++) {
        radi[i] = 0.;
      }

      for (int i = 0; i < res.x * res.y; i++) {
        RGBSpectrum rgb = map[i];
        data[i] = rgb[0];

        float theta = (float) (i / res.x) / res.y * 180.;
        float phi   = (float) (i % res.x) / res.x * 360.;
        radi[ (int) std::sqrt( theta * theta + phi * phi) ] += rgb[0];
      }
      distr1 = new Distribution1D(radi, res.x * res.y);
      distr2 = new Distribution2D(data, res.x, res.y);
    }
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;
  private:
    Distribution1D* distr1;
    Distribution2D* distr2;
};

IceMaterial *CreateIceMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_ICE_H
