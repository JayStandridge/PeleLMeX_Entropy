#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

#define DUMPDATA true


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
smoothed_func( const amrex::Real x[AMREX_SPACEDIM],
	       amrex::Real s_ext[NVAR],
	       ProbParm const& prob_parm,
	       bool isInit)
{
   auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
   

    amrex::Real massfrac_coflow[NUM_SPECIES] = {0.0};

    massfrac_coflow[O2_ID]  = .21;
    massfrac_coflow[N2_ID]  = .78;
    massfrac_coflow[AR_ID]  = .01;

    amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};

    massfrac_jet[CH4_ID] = .1558;
    massfrac_jet[O2_ID]  = .1958;
    massfrac_jet[N2_ID]  = .6367;
    massfrac_jet[AR_ID]  = .0117;
    
    amrex::Real massfrac_pilot[NUM_SPECIES] = {0.0};

    massfrac_pilot[N2_ID]  = .73420;
    massfrac_pilot[O2_ID]  = .05400;
    massfrac_pilot[O_ID]   = .00075;
    massfrac_pilot[H2_ID]  = .00013;
    massfrac_pilot[H_ID]   = .00002;
    massfrac_pilot[H2O_ID] = .09420;
    massfrac_pilot[CO_ID]  = .00410;
    massfrac_pilot[CO2_ID] = .10980;
    massfrac_pilot[OH_ID]  = .00280;
    
   
    
    if (isInit){
      s_ext[TEMP] = prob_parm.T_coflow;
      AMREX_D_TERM(s_ext[VELX] =0.0;,
		   s_ext[VELY] =  prob_parm.U_coflow;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] = massfrac_coflow[n];
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
      return;
    }

    
    // JET //
    if ( x[0] < .003725-prob_parm.R_smooth ) {
      s_ext[TEMP] = prob_parm.T_jet;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = prob_parm.U_jet;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] = massfrac_jet[n];
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }




      // JET to PILOT //
    } else if ( std::abs(x[0]-.003725)<=prob_parm.R_smooth ) {
      amrex::Real val=std::erf(2.5*(x[0]-.003725)/prob_parm.R_smooth)/std::erf(2.5);
      s_ext[TEMP] = (prob_parm.T_jet+prob_parm.T_pilot)/2 -(prob_parm.T_jet-prob_parm.T_pilot)/2*val;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = (prob_parm.U_jet+prob_parm.U_pilot)/2 -(prob_parm.U_jet-prob_parm.U_pilot)/2*val;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] =(massfrac_jet[n]+massfrac_pilot[n])/2 -(massfrac_jet[n]-massfrac_pilot[n])/2*val;
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }



      // PILOT //
    }else if (  x[0] < .009275-prob_parm.R_smooth_2 ) {
      s_ext[TEMP] = prob_parm.T_pilot;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = prob_parm.U_pilot;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] = massfrac_pilot[n];
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }





      // PILOT to COFLOW 1 //
    }else if ( x[0]<=.009275 ) {
      amrex::Real val=std::erf(2.5*(x[0]-.009275)/prob_parm.R_smooth_2)/std::erf(2.5);
      s_ext[TEMP] = (prob_parm.T_pilot+prob_parm.T_coflow)/2 -(prob_parm.T_pilot-prob_parm.T_coflow)/2*val;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = (prob_parm.U_pilot+prob_parm.U_coflow)/2 -(prob_parm.U_pilot-prob_parm.U_coflow)/2*val;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] =(massfrac_pilot[n]+massfrac_coflow[n])/2 -(massfrac_pilot[n]-massfrac_coflow[n])/2*val;
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }


      // JET TO COFLOW 2
    }else if ( x[0]<=.009275+prob_parm.R_smooth_3 ) {
      amrex::Real val=std::erf(2.5*(x[0]-.009275)/prob_parm.R_smooth_3)/std::erf(2.5);
      s_ext[TEMP] = (prob_parm.T_pilot+prob_parm.T_coflow)/2 -(prob_parm.T_pilot-prob_parm.T_coflow)/2*val;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = (prob_parm.U_pilot+prob_parm.U_coflow)/2 -(prob_parm.U_pilot-prob_parm.U_coflow)/2*val;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] =(massfrac_pilot[n]+massfrac_coflow[n])/2 -(massfrac_pilot[n]-massfrac_coflow[n])/2*val;
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }


      // COFLOW
    }else{
      s_ext[TEMP] = prob_parm.T_coflow;
      AMREX_D_TERM(s_ext[VELX] = 0.0;,
		   s_ext[VELY] = prob_parm.U_coflow;,
		   s_ext[VELZ] = 0.0);
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      for (int n = 0; n < NUM_SPECIES; n++) {
	massfrac[n] = massfrac_coflow[n];
      }

      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
	s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
      
      
      
      
    }
}
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real x[AMREX_SPACEDIM] = {0.0};
    AMREX_D_TERM( x[0] = prob_lo[0] + (i+0.5)*dx[0];,
                  x[1] = prob_lo[1] + (j+0.5)*dx[1];,
                  x[2] = prob_lo[2] + (k+0.5)*dx[2];);

    amrex::Real s_ext[NVAR] = {0.0};
    
    smoothed_func(x, s_ext, prob_parm, false);
    
    // amrex::Real xp[AMREX_SPACEDIM] = {0.0};
    // amrex::Real xm[AMREX_SPACEDIM] = {0.0};
    // amrex::Real s_extp[NVAR] = {0.0};
    // amrex::Real s_extm[NVAR] = {0.0};
    // xp[1]=(2*std::pow(x[0],2)*x[1]+.15)/3;
    // xm[1]=2*x[1]/3;
    // smoothed_func(xp, s_extp, prob_parm);
    // smoothed_func(xm, s_extm, prob_parm);
    // amrex::Real val1 = 0.0;
    // amrex::Real val2 = 0.0;
    // val1 = 1-std::pow(x[0],2);
    // val2 = (1-val1)/2;
    for (int n=0; n<NVAR; n++){
      state(i,j,k,n)=s_ext[n];
    }
}

    


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
  smoothed_func(x, s_ext, prob_parm, false);
}



AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);

    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}



#endif
