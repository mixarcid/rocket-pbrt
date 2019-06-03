
// materials/glass.cpp*
#include "materials/ice.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"
#include "interaction.h"

namespace pbrt {

// IceMaterial Method Definitions
void IceMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
					     MemoryArena &arena,
					     TransportMode mode,
					     bool allowMultipleLobes) const {
  si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
  si->bsdf->Add(ARENA_ALLOC(arena, IceBSDF)(distr));
}

IceMaterial *CreateIceMaterial(const TextureParams &mp) {
    std::string mapName = mp.FindFilename("mapname", "");
    return new IceMaterial(mapName);
}

}  // namespace pbrt
