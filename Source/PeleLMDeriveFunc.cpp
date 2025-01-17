#include "PeleLM_Index.H"
#include "PeleLM.H"
#include "PeleLM_K.H"
#include "PeleLMDeriveFunc.H"
#include "pelelm_prob.H"

#include <PelePhysics.H>
#include <mechanism.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

//
// Extract temp
//
void pelelm_dertemp (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geom*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,TEMP);
    });
}

//
// Compute heat release
//
void pelelm_derheatrelease (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                            const FArrayBox& statefab, const FArrayBox& reactfab, const FArrayBox& /*pressfab*/,
                            const Geometry& /*geom*/,
                            Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(!a_pelelm->m_incompressible);

    FArrayBox EnthFab;
    EnthFab.resize(bx,NUM_SPECIES,The_Async_Arena());

    auto const temp = statefab.const_array(TEMP);
    auto const react = reactfab.const_array(0);
    auto const& Hi    = EnthFab.array();
    auto       HRR = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getHGivenT( i, j, k, temp, Hi );
        HRR(i,j,k) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
           HRR(i,j,k) -= Hi(i,j,k,n) * react(i,j,k,n);
        }
    });
}

//
// Extract species mass fractions Y_n
//
void pelelm_dermassfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geom*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,DENSITY);
        der(i,j,k,n) = in_dat(i,j,k,FIRSTSPEC+n) * rhoinv;
    });
}

//
// Extract species mole fractions X_n
//
void pelelm_dermolefrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geom*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real Xt[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,DENSITY);
        for (int n = 0; n < NUM_SPECIES; n++) {
           Yt[n] = in_dat(i,j,k,FIRSTSPEC+n) * rhoinv;
        }
        auto eos = pele::physics::PhysicsType::eos();
        eos.Y2X(Yt,Xt);
        for (int n = 0; n < NUM_SPECIES; n++) {
           der(i,j,k,n) = Xt[n];
        }
    });
}

//
// Extract rho - sum rhoY
//
void pelelm_derrhomrhoy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geom*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = in_dat(i,j,k,DENSITY);
        for (int n = 0; n < NUM_SPECIES; n++) {
            der(i,j,k,0) -= in_dat(i,j,k,FIRSTSPEC+n);
        }
    });
}

//
// Compute cell averaged pressure from nodes
//
void pelelm_deravgpress (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& /*statefab*/, const FArrayBox& /*reactfab*/, const FArrayBox& pressfab,
                         const Geometry& /*geom*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    auto const in_dat = pressfab.array();
    auto       der = derfab.array(dcomp);
    Real factor = 1.0 / ( AMREX_D_TERM(2.0,*2.0,*2.0) );
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)
#if (AMREX_SPACEDIM >= 2 )
                                + in_dat(i+1,j+1,k)   + in_dat(i,j+1,k)
#if (AMREX_SPACEDIM == 3 )
                                + in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)
                                + in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1)
#endif
#endif
                                );
    });
}

//
// Compute the velocity magnitude
//
void pelelm_dermgvel (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                      const Geometry& /*geom*/,
                      Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    auto const vel = statefab.array(VELX);
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = std::sqrt(( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                              + vel(i,j,k,1)*vel(i,j,k,1),
                                              + vel(i,j,k,2)*vel(i,j,k,2)) ));
    });
}

//
// Compute vorticity magnitude
//
void pelelm_dermgvort (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                       const Geometry& geom,
                       Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);,
                 const amrex::Real idy = geom.InvCellSize(1);,
                 const amrex::Real idz = geom.InvCellSize(2););

    auto const& dat_arr = statefab.const_array();
    auto const&vort_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(statefab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

    auto typ = flags.getType(bx);

    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            vort_arr(i,j,k) = 0.0;
        });
    } else if (typ == FabType::singlevalued) {
        const auto& flag_fab = flags.const_array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real c0 = -1.5;
            constexpr amrex::Real c1 = 2.0;
            constexpr amrex::Real c2 = -0.5;
            if (flag_fab(i,j,k).isCovered()) {
               vort_arr(i,j,k) = 0.0;
            } else {
               // Define interpolation lambda
               auto onesided = [] (const Real &v0,
                                   const Real &v1,
                                   const Real &v2) -> Real
               {return c0*v0+c1*v1+c2*v2;};

               amrex::Real vx = 0.0;
               amrex::Real uy = 0.0;
#if ( AMREX_SPACEDIM == 2 )
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
               } else {
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected( 0,1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
               }
               vort_arr(i,j,k) = std::abs(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
               amrex::Real wx = 0.0;
               amrex::Real wy = 0.0;
               amrex::Real uz = 0.0;
               amrex::Real vz = 0.0;
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   // Covered cell to the right, go fish left
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                   wx = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   // Covered cell to the left, go fish right
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                   wx = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
               } else {
                   // No covered cells right or left, use standard stencil
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                   wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                   wy = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                   wy = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                   wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
               }
               // Do the same in z-direction
               if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                   uz = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                   vz = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
               } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                   uz = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                   vz = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
               } else {
                   uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                   vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
               }
               vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
            }
        });
    } else
#endif          // Check on EB
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            vort_arr(i,j,k) = std::abs(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }
}

//
// Compute vorticity components
//
void pelelm_dervort (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& geom,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);,
                 const amrex::Real idy = geom.InvCellSize(1);,
                 const amrex::Real idz = geom.InvCellSize(2););

    auto const& dat_arr = statefab.const_array();
    auto const&vort_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(statefab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

    auto typ = flags.getType(bx);

    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            vort_arr(i,j,k,n) = 0.0;
        });
    } else if (typ == FabType::singlevalued) {
        const auto& flag_fab = flags.const_array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real c0 = -1.5;
            constexpr amrex::Real c1 = 2.0;
            constexpr amrex::Real c2 = -0.5;
            if (flag_fab(i,j,k).isCovered()) {
                for (int n{0}; n < ncomp; ++n) {
                    vort_arr(i,j,k,n) = 0.0;
                }
            } else {
               // Define interpolation lambda
               auto onesided = [] (const Real &v0,
                                   const Real &v1,
                                   const Real &v2) -> Real
               {return c0*v0+c1*v1+c2*v2;};

               amrex::Real vx = 0.0;
               amrex::Real uy = 0.0;
#if ( AMREX_SPACEDIM == 2 )
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
               } else {
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected( 0,1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
               }
               vort_arr(i,j,k) = vx-uy;

#elif ( AMREX_SPACEDIM == 3 )
               amrex::Real wx = 0.0;
               amrex::Real wy = 0.0;
               amrex::Real uz = 0.0;
               amrex::Real vz = 0.0;
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   // Covered cell to the right, go fish left
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                   wx = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   // Covered cell to the left, go fish right
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                   wx = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
               } else {
                   // No covered cells right or left, use standard stencil
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                   wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                   wy = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                   wy = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                   wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
               }
               // Do the same in z-direction
               if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                   uz = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                   vz = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
               } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                   uz = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                   vz = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
               } else {
                   uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                   vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
               }
               vort_arr(i,j,k,0) = (wy-vz)*(wy-vz);
               vort_arr(i,j,k,1) = (uz-wx)*(uz-wx);
               vort_arr(i,j,k,2) = (vx-uy)*(vx-uy);
#endif
            }
        });
    } else
#endif          // Check on EB

    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            vort_arr(i,j,k) = vx-uy;

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k,0) = (wy-vz)*(wy-vz);
            vort_arr(i,j,k,1) = (uz-wx)*(uz-wx);
            vort_arr(i,j,k,2) = (vx-uy)*(vx-uy);
#endif
        });
    }
}

//
// Compute cell-centered coordinates
//
void pelelm_dercoord (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                      const Geometry& geom,
                      Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_D_TERM(const amrex::Real dx = geom.CellSize(0);,
                 const amrex::Real dy = geom.CellSize(1);,
                 const amrex::Real dz = geom.CellSize(2););

    auto const& coord_arr = derfab.array(dcomp);
    const auto geomdata = geom.data();

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(statefab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

    auto typ = flags.getType(bx);
    // Compute cell center coordinates even in covered boxes/cell. Only
    // modify the cell-center in cut cells
    if (typ == FabType::singlevalued) {
        const auto& flag_arr = flags.const_array();
        const auto& ccent_fab = ebfab.getCentroidData();
        const auto& ccent_arr = ccent_fab->const_array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real* prob_lo = geomdata.ProbLo();
            if (flag_arr(i,j,k).isCovered() || flag_arr(i,j,k).isRegular()) {
                AMREX_D_TERM(coord_arr(i,j,k,0) = prob_lo[0] + (i+0.5)*dx;,
                             coord_arr(i,j,k,1) = prob_lo[1] + (j+0.5)*dy;,
                             coord_arr(i,j,k,2) = prob_lo[2] + (k+0.5)*dz;);
            } else {
                AMREX_D_TERM(coord_arr(i,j,k,0) = prob_lo[0] + (i+0.5+ccent_arr(i,j,k,0))*dx;,
                             coord_arr(i,j,k,1) = prob_lo[1] + (j+0.5+ccent_arr(i,j,k,1))*dy;,
                             coord_arr(i,j,k,2) = prob_lo[2] + (k+0.5+ccent_arr(i,j,k,2))*dz;);
            }
        });
    } else
#endif
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real* prob_lo = geomdata.ProbLo();
            AMREX_D_TERM(coord_arr(i,j,k,0) = prob_lo[0] + (i+0.5)*dx;,
                         coord_arr(i,j,k,1) = prob_lo[1] + (j+0.5)*dy;,
                         coord_arr(i,j,k,2) = prob_lo[2] + (k+0.5)*dz;);
        });
    }
}

//
// Compute Q-criterion
//
void pelelm_derQcrit (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                      const Geometry& geom,
                      Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
#if AMREX_SPACEDIM == 3
    AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);,
                 const amrex::Real idy = geom.InvCellSize(1);,
                 const amrex::Real idz = geom.InvCellSize(2););

    auto const &  dat_arr = statefab.const_array();
    auto const &qcrit_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(statefab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

    auto typ = flags.getType(bx);

    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            qcrit_arr(i,j,k) = 0.0;
        });
    } else if (typ == FabType::singlevalued) {
        const auto& flag_fab = flags.const_array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag_fab(i,j,k).isCovered()) {
                qcrit_arr(i,j,k) = 0.0;
            } else {
                // Define interpolation lambda
                constexpr amrex::Real c0 = -1.5;
                constexpr amrex::Real c1 = 2.0;
                constexpr amrex::Real c2 = -0.5;
                auto onesided = [] (const Real &v0,
                                    const Real &v1,
                                    const Real &v2) -> Real
                {return c0*v0+c1*v1+c2*v2;};

                // Strain rate tensor
                Array2D<Real,0,2,0,2> gradU;
                if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                    gradU(0,0) = - onesided(dat_arr(i,j,k,0), dat_arr(i-1,j,k,0), dat_arr(i-2,j,k,0)) * idx;
                    gradU(1,0) = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                    gradU(2,0) = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
                } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                    gradU(0,0) = onesided(dat_arr(i,j,k,0), dat_arr(i+1,j,k,0), dat_arr(i+2,j,k,0)) * idx;
                    gradU(1,0) = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                    gradU(2,0) = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
                } else {
                    gradU(0,0) = 0.5 * (dat_arr(i+1,j,k,0) - dat_arr(i-1,j,k,0)) * idx;
                    gradU(1,0) = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                    gradU(2,0) = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
                }
                if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                    gradU(0,1) = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                    gradU(1,1) = - onesided(dat_arr(i,j,k,1), dat_arr(i,j-1,k,1), dat_arr(i,j-2,k,1)) * idy;
                    gradU(2,1) = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
                } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                    gradU(0,1) = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                    gradU(1,1) = onesided(dat_arr(i,j,k,1), dat_arr(i,j+1,k,1), dat_arr(i,j+2,k,1)) * idy;
                    gradU(2,1) = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
                } else {
                    gradU(0,1) = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                    gradU(1,1) = 0.5 * (dat_arr(i,j+1,k,1) - dat_arr(i,j-1,k,1)) * idy;
                    gradU(2,1) = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
                }
                if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                    gradU(0,2) = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                    gradU(1,2) = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
                    gradU(2,2) = - onesided(dat_arr(i,j,k,2), dat_arr(i,j,k-1,2), dat_arr(i,j,k-2,2)) * idz;
                } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                    gradU(0,2) = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                    gradU(1,2) = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
                    gradU(2,2) = onesided(dat_arr(i,j,k,2), dat_arr(i,j,k+1,2), dat_arr(i,j,k+2,2)) * idz;
                } else {
                    gradU(0,2) = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                    gradU(1,2) = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
                    gradU(2,2) = 0.5 * (dat_arr(i,j,k+1,2) - dat_arr(i,j,k-1,2)) * idz;
                }

                // Divu
                amrex::Real divU = gradU(0,0) + gradU(1,1) + gradU(2,2);

                // Directly Assemble Sym. & AntiSym. into Qcrit.
                // Remove divU (dilatation) from the Sym. tensor (due to mixing/reaction most often)
                qcrit_arr(i,j,k) = 0.0;
                for (int dim1 = 0; dim1 < AMREX_SPACEDIM; ++dim1) {
                  for (int dim2 = 0; dim2 < AMREX_SPACEDIM; ++dim2) {
                    Real Ohm = 0.5 * (gradU(dim1,dim2) - gradU(dim2,dim1));
                    Real Sij = 0.5 * (gradU(dim1,dim2) + gradU(dim2,dim1));
                    if (dim1 == dim2) {
                      Sij -= divU/AMREX_SPACEDIM;
                    }
                    qcrit_arr(i,j,k) += Ohm*Ohm - Sij*Sij;
                  }
                }
            }
        });
    } else
#endif
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // Strain rate tensor
          Array2D<Real,0,2,0,2> gradU;
          gradU(0,0) = 0.5 * (dat_arr(i+1,j,k,0) - dat_arr(i-1,j,k,0)) * idx;
          gradU(0,1) = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
          gradU(0,2) = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
          gradU(1,0) = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
          gradU(1,1) = 0.5 * (dat_arr(i,j+1,k,1) - dat_arr(i,j-1,k,1)) * idy;
          gradU(1,2) = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
          gradU(2,0) = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
          gradU(2,1) = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
          gradU(2,2) = 0.5 * (dat_arr(i,j,k+1,2) - dat_arr(i,j,k-1,2)) * idz;

          // Divu
          amrex::Real divU = gradU(0,0) + gradU(1,1) + gradU(2,2);

          // Directly Assemble Sym. & AntiSym. into Qcrit.
          // Remove divU (dilatation) from the Sym. tensor (due to mixing/reaction most often)
          qcrit_arr(i,j,k) = 0.0;
          for (int dim1 = 0; dim1 < AMREX_SPACEDIM; ++dim1) {
            for (int dim2 = 0; dim2 < AMREX_SPACEDIM; ++dim2) {
              Real Ohm = 0.5 * (gradU(dim1,dim2) - gradU(dim2,dim1));
              Real Sij = 0.5 * (gradU(dim1,dim2) + gradU(dim2,dim1));
              if (dim1 == dim2) {
                Sij -= divU/AMREX_SPACEDIM;
              }
              qcrit_arr(i,j,k) += Ohm*Ohm - Sij*Sij;
            }
          }
        });
    }
#endif

}

//
// Compute the kinetic energy
//
void pelelm_derkineticenergy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                              const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                              const Geometry& /*geom*/,
                              Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    if (a_pelelm->m_incompressible) {
       auto const vel = statefab.array(VELX);
       auto       der = derfab.array(dcomp);
       amrex::ParallelFor(bx,
       [=,rho=a_pelelm->m_rho] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           der(i,j,k) = 0.5 * rho
                            * ( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                             + vel(i,j,k,1)*vel(i,j,k,1),
                                             + vel(i,j,k,2)*vel(i,j,k,2)) );
       });
    } else {
       auto const rho = statefab.array(DENSITY);
       auto const vel = statefab.array(VELX);
       auto       der = derfab.array(dcomp);
       amrex::ParallelFor(bx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           der(i,j,k) = 0.5 * rho(i,j,k)
                            * ( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                             + vel(i,j,k,1)*vel(i,j,k,1),
                                             + vel(i,j,k,2)*vel(i,j,k,2)) );
       });
    }
}

//
// Compute enstrophy
//
void pelelm_derenstrophy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                          const Geometry& geom,
                          Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);,
                 const amrex::Real idy = geom.InvCellSize(1);,
                 const amrex::Real idz = geom.InvCellSize(2););

    auto const&  dat_arr = statefab.const_array(VELX);
    auto const&  rho_arr = (a_pelelm->m_incompressible) ? Array4<const Real>{}
                                                        : statefab.const_array(DENSITY);
    auto const&  ens_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
    const EBFArrayBox& ebfab = static_cast<EBFArrayBox const&>(statefab);
    const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

    auto typ = flags.getType(bx);

    if (typ == FabType::covered)
    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ens_arr(i,j,k) = 0.0;
        });
    } else if (typ == FabType::singlevalued) {
        const auto& flag_fab = flags.const_array();
        amrex::ParallelFor(bx, [=,incomp=a_pelelm->m_incompressible,rho=a_pelelm->m_rho]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real c0 = -1.5;
            constexpr amrex::Real c1 = 2.0;
            constexpr amrex::Real c2 = -0.5;
            if (flag_fab(i,j,k).isCovered()) {
                ens_arr(i,j,k) = 0.0;
            } else {
               Real l_rho = rho;
               if (!incomp) {
                   l_rho = rho_arr(i,j,k);
               }
               // Define interpolation lambda
               auto onesided = [] (const Real &v0,
                                   const Real &v1,
                                   const Real &v2) -> Real
               {return c0*v0+c1*v1+c2*v2;};

               amrex::Real vx = 0.0;
               amrex::Real uy = 0.0;
#if ( AMREX_SPACEDIM == 2 )
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
               } else {
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected( 0,1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
               }
               ens_arr(i,j,k) = 0.5 * l_rho * (vx-uy)*(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
               amrex::Real wx = 0.0;
               amrex::Real wy = 0.0;
               amrex::Real uz = 0.0;
               amrex::Real vz = 0.0;
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   // Covered cell to the right, go fish left
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                   wx = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   // Covered cell to the left, go fish right
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                   wx = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
               } else {
                   // No covered cells right or left, use standard stencil
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                   wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                   wy = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                   wy = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                   wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
               }
               // Do the same in z-direction
               if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                   uz = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                   vz = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
               } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                   uz = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                   vz = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
               } else {
                   uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                   vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
               }
               ens_arr(i,j,k) = 0.5 * l_rho * ((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
            }
        });
    } else
#endif
    {
        amrex::ParallelFor(bx, [=,incomp=a_pelelm->m_incompressible,rho=a_pelelm->m_rho]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real l_rho = rho;
            if (!incomp) {
                l_rho = rho_arr(i,j,k);
            }
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            ens_arr(i,j,k) = 0.5 * l_rho * (vx-uy)*(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            ens_arr(i,j,k) = 0.5 * l_rho * ((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }
}


//
// Compute mixture fraction
//
void pelelm_dermixfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geom*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1);

    if (a_pelelm->Zfu < 0.0) amrex::Abort("Mixture fraction not initialized");

    auto const density   = statefab.array(DENSITY);
    auto const rhoY      = statefab.array(FIRSTSPEC);
    auto       mixt_frac = derfab.array(dcomp);

    amrex::Real Zox_lcl   = a_pelelm->Zox;
    amrex::Real Zfu_lcl   = a_pelelm->Zfu;
    amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_Bilger;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_Bilger[n] = a_pelelm->spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, rhoY, mixt_frac, fact_Bilger, Zox_lcl, denom_inv] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        mixt_frac(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) += ( rhoY(i,j,k,n) * fact_Bilger[n] ) * rho_inv;
        }
        mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - Zox_lcl ) * denom_inv;
    });
}

//
// Compute progress variable
//
void pelelm_derprogvar (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geom*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1);

    if (a_pelelm->m_C0 < 0.0) amrex::Abort("Progress variable not initialized");

    auto const density  = statefab.array(DENSITY);
    auto const rhoY     = statefab.array(FIRSTSPEC);
    auto const temp     = statefab.array(TEMP);
    auto       prog_var = derfab.array(dcomp);

    amrex::Real C0_lcl   = a_pelelm->m_C0;
    amrex::Real C1_lcl   = a_pelelm->m_C1;
    amrex::Real denom_inv = 1.0 / (C1_lcl - C0_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES+1> Cweights;
    for (int n=0; n<NUM_SPECIES+1; ++n) {
        Cweights[n] = a_pelelm->m_Cweights[n];
    }

    amrex::ParallelFor(bx, [=,revert=a_pelelm->m_Crevert]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        prog_var(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            prog_var(i,j,k) += ( rhoY(i,j,k,n) * Cweights[n] ) * rho_inv;
        }
        prog_var(i,j,k) += temp(i,j,k) * Cweights[NUM_SPECIES];
        if (revert) {
           prog_var(i,j,k) = 1.0 - ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        } else {
           prog_var(i,j,k) = ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        }
    });
}

//
// Extract mixture viscosity
//
void pelelm_dervisc (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geom*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

    if (a_pelelm->m_incompressible) {
        derfab.setVal<RunOn::Device>(a_pelelm->m_mu,bx,dcomp,1);
    } else {
        auto const& rhoY = statefab.const_array(FIRSTSPEC);
        auto const& T    = statefab.array(TEMP);
        auto       der   = derfab.array(dcomp);
        auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
        amrex::ParallelFor(bx,
        [rhoY,T,der,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            getVelViscosity(i, j, k, rhoY, T, der, ltransparm);
        });
    }
}

//
// Extract mixture averaged species diffusion coefficients
//
void pelelm_derdiffc (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                      const Geometry& /*geom*/,
                      Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    if (a_pelelm->m_use_soret) {
      AMREX_ASSERT(ncomp == 2*NUM_SPECIES);
    } else {
      AMREX_ASSERT(ncomp == NUM_SPECIES);
    }
    FArrayBox dummies(bx,2,The_Async_Arena());
    auto const& rhoY = statefab.const_array(FIRSTSPEC);
    auto const& T    = statefab.array(TEMP);
    auto       rhoD  = derfab.array(dcomp);
    auto     lambda  = dummies.array(0);
    auto         mu  = dummies.array(1);
    auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
    if (a_pelelm->m_use_soret) {
      auto rhotheta = derfab.array(dcomp+NUM_SPECIES);
      amrex::ParallelFor(bx,
        [rhoY,T,rhoD,rhotheta,lambda,mu,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          getTransportCoeffSoret(i, j, k, rhoY, T, rhoD, rhotheta, lambda, mu, ltransparm);
        });
    } else {
      if (a_pelelm->m_unity_Le) {
        amrex::Real ScInv = a_pelelm->m_Schmidt_inv;
        amrex::Real PrInv = a_pelelm->m_Prandtl_inv;
        amrex::ParallelFor(bx,
          [rhoY,T,rhoD,lambda,mu,ltransparm,ScInv,PrInv] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            getTransportCoeffUnityLe(i, j, k, ScInv, PrInv, rhoY, T, rhoD, lambda, mu, ltransparm);
          });
      } else {
        amrex::ParallelFor(bx,
          [rhoY,T,rhoD,lambda,mu,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            getTransportCoeff(i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
          });
      }
    }
}

//
// Extract thermal diffusivity
//
void pelelm_derlambda (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                       const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                       const Geometry& /*geom*/,
                       Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

    FArrayBox dummies(bx,NUM_SPECIES+1,The_Async_Arena());
    auto const& rhoY = statefab.const_array(FIRSTSPEC);
    auto const& T    = statefab.array(TEMP);
    auto       rhoD  = dummies.array(1);
    auto     lambda  = derfab.array(dcomp);
    auto         mu  = dummies.array(0);
    auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
    amrex::Real ScInv = a_pelelm->m_Schmidt_inv;
    amrex::Real PrInv = a_pelelm->m_Prandtl_inv;
    int unity_Le = a_pelelm->m_unity_Le;
    amrex::ParallelFor(bx,
    [rhoY,T,rhoD,lambda,mu,ltransparm,unity_Le,ScInv,PrInv]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      if (unity_Le) {
        getTransportCoeffUnityLe(i, j, k, ScInv, PrInv, rhoY, T, rhoD, lambda, mu, ltransparm);
      } else {
        getTransportCoeff(i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
      }
    });
}

//
// Extract Distribution Mapping
//
void pelelm_derdmap (PeleLM* /*a_pelelm*/, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& /*statefab*/, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geom*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    auto       der = derfab.array(dcomp);
    const int myrank = ParallelDescriptor::MyProc();
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = myrank;
    });
}

//
// Compute EITERM4
//
void pelelm_derEIterm4 (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                            const FArrayBox& statefab, const FArrayBox& reactfab, const FArrayBox& /*pressfab*/,
                            const Geometry& geomdata,
                            Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(!a_pelelm->m_incompressible);

    FArrayBox EnthFab;
    EnthFab.resize(bx,NUM_SPECIES,The_Async_Arena());

    auto const dat = statefab.const_array(); 
    auto const react = reactfab.const_array(0);
    auto const& Hi    = EnthFab.array();
    auto       EI = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
      EI(i,j,k) = 0.0;
      amrex::Real massfrac[NUM_SPECIES];
      amrex::Real gibbs_fe[NUM_SPECIES] = {0.0};
      amrex::Real gibbs_standard[NUM_SPECIES] = {0.0};
      amrex::Real prod_rate[NUM_SPECIES] = {0.0};
      
      amrex::Real rho = dat(i, j, k, DENSITY);
      amrex::Real rhoInv = 1/rho;
      // Get rate of production of each species
      amrex::Real sc[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] = dat(i, j, k, FIRSTSPEC + n) * rhoInv;
        sc[n] = rho * massfrac[n] *imw(n);
      }
      CKYTCR(rho,rho,massfrac, sc); //Need species concentration first
      
      // Get chemical potential, or molar gibbs
      
      amrex::Real tc[5]={0.0};
      tc[0]= std::log(dat(i,j,k,TEMP));
      tc[1]=dat(i,j,k,TEMP);
      tc[2]=tc[1]*tc[1];
      tc[3]=tc[2]*tc[1];
      tc[4]=tc[3]*tc[1];
      gibbs(gibbs_fe, tc);
      for (int n = 0; n < NUM_SPECIES; n++) {
	gibbs_standard[n] = 8.31446 * dat(i,j,k,TEMP) * gibbs_fe[n];
	gibbs_fe[n] += std::log(std::max(sc[n]*8.31446 * dat(i,j,k,TEMP)/101325,std::numeric_limits<double>::epsilon()));
	gibbs_fe[n] *= 8.31446 * dat(i,j,k,TEMP);
	
      }
      for (int ii = 0; ii < NUM_SPECIES; ii++) { 
	sc[ii] *= 1e6; // in SI units for productionRate
      }
      productionRate(prod_rate, sc, dat(i,j,k,TEMP));
      EI(i,j,k) = 0.0;
      
      // for (int ii = 0; ii < NUM_SPECIES; ii++) { 
      // 	EI(i,j,k)+=prod_rate[ii]*gibbs_fe[ii];
	
      // }
      amrex::Real q_f_temp[NUM_REACTIONS] = {0.0};
      amrex::Real q_r_temp[NUM_REACTIONS] = {0.0};
      amrex::Real q_f[NUM_REACTIONS] = {0.0};
      amrex::Real q_r[NUM_REACTIONS] = {0.0};
      amrex::Real prod_rate_2[NUM_SPECIES] = {0.0};
      amrex::Real wdot[NUM_REACTIONS] = {0.0};
      amrex::Real nu_spec = 0.0;
      int nspec = 0;
      int& nspecr = nspec;
      int* temp;
      int* temp2;
      CKINU(0, nspecr, temp, temp2);
      int ki[nspec] = {0};
      int nu[nspec] = {0};
      int* kip = ki;
      int* nup = nu;
      int rmap[NUM_REACTIONS] = {0};
      int* rmapp = rmap;
      GET_RMAP(rmapp);
      for (int ii = 0; ii < NUM_SPECIES; ii++) { 
      	sc[ii] *= 1e-6; // in SI units for productionRate
      }
      comp_qfqr(q_f_temp, q_r_temp, sc,sc,tc,1.0/tc[1]);
      
      for (int n = 0; n < NUM_REACTIONS; n++) { 
       
       	q_f[rmap[n]] = q_f_temp[n];
       	q_r[rmap[n]] = q_r_temp[n];

      }
      double DG_j[NUM_REACTIONS] = {0.0};
      double DG_j_s[NUM_REACTIONS] = {0.0};
      double EI_j[NUM_REACTIONS] = {0.0};
      bool skip = false;
      double EImin = 0.0;
      int EIminInd = 0;
      for (int n = 0; n < NUM_REACTIONS; n++) {
	skip = false;
      	CKINU(n+1, nspecr, kip, nup);
      	wdot[n] = 1e-6 * (q_f[n] - q_r[n]);
      	for (int m = 0; m < nspec; m++){
      	  DG_j[n] += nu[m] * gibbs_fe[ki[m]-1];
      	  DG_j_s[n] += nu[m] * gibbs_standard[ki[m]-1];
	  if (sc[ki[m]-1] < 0.0){
	    skip = true;
	  }
      	}
	
	EI_j[n] = wdot[n] * DG_j[n];
	
	if (!skip){
	  if (EI_j[n] < EImin){
	    EImin = EI_j[n];
	    EIminInd = n;
	  }
	  EI(i,j,k) += EI_j[n];
	}
      }
#if DUMPDATA==true
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx      = geomdata.CellSize();
      
      auto eos = pele::physics::PhysicsType::eos();
      amrex::Real x[AMREX_SPACEDIM] = {0.0};
      AMREX_D_TERM( x[0] = prob_lo[0] + (i+0.5)*dx[0];,
		    x[1] = prob_lo[1] + (j+0.5)*dx[1];,
		    x[2] = prob_lo[2] + (k+0.5)*dx[2];);
      
      // Open the file in append mode
      std::ofstream file("data.txt", std::ios::app);
      
      // Check if the file is opened successfully
      if (!file.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
      }
      
      // Append the data row to the file
      file << x[0] << ", " << x[1] << ", " << dat(i,j,k,DENSITY) << ", " << dat(i,j,k,TEMP)<< std::endl;
      
      // Close the file
      file.close();

#endif
    });
}

