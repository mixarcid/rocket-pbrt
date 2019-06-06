// media/grid.cpp*
#include "media/emissive.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

namespace pbrt {

// EmissiveMedium Method Definitions
Float EmissiveMedium::Density(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // Trilinearly interpolate density values to compute local density
    Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
    Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
    Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
    Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
    Float d0 = Lerp(d.y, d00, d10);
    Float d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

Spectrum EmissiveMedium::getLe(const Point3f &p) const {
    // Compute voxel coordinates and offsets for _p_
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    Spectrum le1 = LeHelper(pi);
    Spectrum le2 = LeHelper(pi + Vector3i(1, 0, 0));
    Spectrum le3 = LeHelper(pi + Vector3i(0, 1, 0));
    Spectrum le4 = LeHelper(pi + Vector3i(1, 1, 0));
    Spectrum le5 = LeHelper(pi + Vector3i(0, 0, 1));
    Spectrum le6 = LeHelper(pi + Vector3i(1, 0, 1));
    Spectrum le7 = LeHelper(pi + Vector3i(0, 1, 1));
    Spectrum le8 = LeHelper(pi + Vector3i(1, 1, 1));

    Spectrum ret;
    for (int i=0; i<3; ++i) {
      Float d00 = Lerp(d.x, le1[i], le2[i]);
      Float d10 = Lerp(d.x, le3[i], le4[i]);
      Float d01 = Lerp(d.x, le5[i], le6[i]);
      Float d11 = Lerp(d.x, le7[i], le8[i]);
      Float d0 = Lerp(d.y, d00, d10);
      Float d1 = Lerp(d.y, d01, d11);
      ret[i] = Lerp(d.z, d0, d1);
    }
    return ret;
}

  
Spectrum EmissiveMedium::Sample(const Ray &rWorld, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // Run delta-tracking iterations to sample a medium interaction
    Float t = tMin;
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        if (Density(ray(t)) * invMaxDensity > sampler.Get1D()) {
            // Populate _mi_ with medium interaction information and return
            PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
            *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                                    phase);
            return sigma_s / sigma_t;
        }
    }
    return Spectrum(1.f);
}

Spectrum EmissiveMedium::Tr(const Ray &rWorld, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);

    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // Perform ratio tracking to estimate the transmittance value
    Float Tr = 1, t = tMin;
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        Float density = Density(ray(t));
        Tr *= 1 - std::max((Float)0, density * invMaxDensity);
        // Added after book publication: when transmittance gets low,
        // start applying Russian roulette to terminate sampling.
        const Float rrThreshold = .1;
        if (Tr < rrThreshold) {
            Float q = std::max((Float).05, 1 - Tr);
            if (sampler.Get1D() < q) return 0;
            Tr /= 1 - q;
        }
    }
    return Spectrum(Tr);
}

const int N_ITERS = 100;
Spectrum EmissiveMedium::SampleLe(const Ray &rWorld, Sampler &sampler,
				  MemoryArena &arena,
				  MediumInteraction *mi) const {
  
  Ray ray = WorldToMedium(Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
  // Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
  const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
  Float tMin, tMax;
  if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(0.f);
  //tMax = tMin + sampler.Get1D()*(tMax - tMin);
    
  Float dt = (tMax-tMin)/float(N_ITERS);
  Float dd = dt*ray.d.Length();
  Spectrum L(0);
  for (int i=0; i<N_ITERS; ++i) {
    Float t = tMin + dt*i;
    Float density = Density(ray(t));
    Float d = std::max((Float)0, density * invMaxDensity);
    L = L - sigma_t*(1-d)*L*dd + (getLe(ray(t))*d)*dd;
  }
  return L;
}

}  // namespace pbrt
