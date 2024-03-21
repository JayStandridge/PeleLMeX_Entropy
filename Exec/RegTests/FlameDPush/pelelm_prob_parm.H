#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Shared params
  amrex::Real T_coflow = 291.0_rt;
  amrex::Real P_mean = 100615.725_rt;
  amrex::Real T_pilot = 1880_rt;
  amrex::Real T_jet = 294_rt;
  amrex::Real U_coflow  = 1.0_rt;
  amrex::Real U_pilot  = 15.93_rt;
  amrex::Real U_jet = 63.1_rt;
  amrex::Real R_smooth = .002_rt;
  amrex::Real R_smooth_2 = .002_rt;
  amrex::Real R_smooth_3 = .002_rt;
};
#endif
