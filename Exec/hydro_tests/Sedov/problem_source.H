#ifndef problem_source_H
#define problem_source_H

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{
    using namespace problem;

    auto coord = geomdata.Coord();
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    // Test for species source: create a variable that mirrors the internal energy
    // update it accordingly with `-p div{U}`
    // (See Castro/Source/sources/Castro_thermo.cpp)
    Real divU = 0;
    // radius for non-Cartesian
    Real rp = prob_lo[0] + (static_cast<Real>(i) + 1.5_rt)*dx[0];
    Real rm = prob_lo[0] + (static_cast<Real>(i) - 0.5_rt)*dx[0];
    Real r = 0.5_rt*(rm + rp);
    // compute div{U}
    if (coord == 0) {
        divU = 0.5_rt*(state(i+1,j,k,UMX)/state(i+1,j,k,URHO)  -
                                state(i-1,j,k,UMX)/state(i-1,j,k,URHO))/dx[0];
    } else if (coord == 1) {
        // axisymmetric
        divU = 0.5_rt*(rp*state(i+1,j,k,UMX)/state(i+1,j,k,URHO) -
                                rm*state(i-1,j,k,UMX)/state(i-1,j,k,URHO))/(r*dx[0]);
    } else if (coord == 2) {
        // spherical
        divU = 0.5_rt*(rp*rp*state(i+1,j,k,UMX)/state(i+1,j,k,URHO) -
                                rm*rm*state(i-1,j,k,UMX)/state(i-1,j,k,URHO))/(r*r*dx[0]);
    }
#if AMREX_SPACEDIM >= 2
    divU += 0.5_rt*(state(i,j+1,k,UMY)/state(i,j+1,k,URHO) -
                                state(i,j-1,k,UMY)/state(i,j-1,k,URHO))/dx[1];
#endif
#if AMREX_SPACEDIM == 3
    divU += 0.5_rt*(state(i,j,k+1,UMZ)/state(i,j,k+1,URHO) -
                                state(i,j,k-1,UMZ)/state(i,j,k-1,URHO))/dx[2];
#endif
    // Compute pressure
    Real rho_e = state(i,j,k,UFX);
    Real p = 2./3 * rho_e;
    src(i,j,k,UFX) = -p * divU;
}

#endif //problem_source_H