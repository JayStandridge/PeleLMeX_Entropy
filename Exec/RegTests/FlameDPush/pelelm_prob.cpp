#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   pp.query("T_coflow", prob_parm->T_coflow);
   pp.query("P_mean", prob_parm->P_mean);
   pp.query("T_pilot", prob_parm->T_pilot);
   pp.query("T_jet", prob_parm->T_jet);
   pp.query("U_coflow", prob_parm->U_coflow);
   pp.query("U_pilot", prob_parm->U_pilot);
   pp.query("U_jet", prob_parm->U_jet);
   pp.query("R_smooth", prob_parm->R_smooth);
   pp.query("R_smooth_2", prob_parm->R_smooth_2);
   pp.query("R_smooth_3", prob_parm->R_smooth_3);
}
