
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MEDIA_EMISSIVE_H
#define PBRT_MEDIA_EMISSIVE_H

// media/grid.h*
#include "medium.h"
#include "transform.h"
#include "stats.h"

namespace pbrt {

// EmissiveMedium Declarations
class EmissiveMedium : public Medium {
  public:
    // EmissiveMedium Public Methods
    EmissiveMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      int nx, int ny, int nz, const Transform &mediumToWorld,
		      const Float *d, const Float* le)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          g(g),
          nx(nx),
          ny(ny),
          nz(nz),
          WorldToMedium(Inverse(mediumToWorld)),
          density(new Float[nx * ny * nz]),
          Le(new Spectrum[nx * ny * nz]) {
        memcpy((Float *) density.get(), d, sizeof(Float) * nx * ny * nz);
	memcpy((Spectrum *) Le.get(), le, sizeof(Spectrum) * nx * ny * nz);
        // Precompute values for Monte Carlo sampling of _EmissiveMedium_
        sigma_t = (sigma_a + sigma_s)[0];
        if (Spectrum(sigma_t) != sigma_a + sigma_s)
            Error(
                "EmissiveMedium requires a spectrally uniform attenuation "
                "coefficient!");
        Float maxDensity = 0;
        for (int i = 0; i < nx * ny * nz; ++i)
            maxDensity = std::max(maxDensity, density[i]);
        invMaxDensity = 1 / maxDensity;
    }

    Float Density(const Point3f &p) const;
    Float D(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return density[(p.z * ny + p.y) * nx + p.x];
    }
    Spectrum getLe(const Point3f &p) const;
    Spectrum LeHelper(const Point3i &p) const {
      Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return Le[(p.z * ny + p.y) * nx + p.x];
    }
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;
    Spectrum SampleLe(const Ray &ray, Sampler &sampler, MemoryArena &arena,
		      MediumInteraction *mi) const;

  private:
    // EmissiveMedium Private Data
    const Spectrum sigma_a, sigma_s;
    const Float g;
    const int nx, ny, nz;
    const Transform WorldToMedium;
    std::unique_ptr<Float[]> density;
    std::unique_ptr<Spectrum[]> Le;
    Float sigma_t;
    Float invMaxDensity;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_EMISSIVE_H
