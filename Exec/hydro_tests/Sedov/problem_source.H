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

    // we now need the pressure -- we will assume that the
    // temperature is consistent with the input state
    eos_rep_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.e = state(i,j,k,UFX) / state(i,j,k,URHO); // internal energy per unit volume
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n)/state(i,j,k,URHO);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state(i,j,k,UFX+n)/state(i,j,k,URHO);
    }
#endif
    eos(eos_input_re, eos_state);

    // final source term, -p div{U}
    src(i,j,k,UFX) = - eos_state.p * divU;
}

#endif //problem_source_H