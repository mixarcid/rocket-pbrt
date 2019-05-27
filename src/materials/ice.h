#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_ICE_H
#define PBRT_MATERIALS_ICE_H

// materials/glass.h*
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// GlassMaterial Declarations
class IceMaterial : public Material {
  public:
    // GlassMaterial Public Methods
    GlassMaterial() {}
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;
};

IceMaterial *CreateIceMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_ICE_H
