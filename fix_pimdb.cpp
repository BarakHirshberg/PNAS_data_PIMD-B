/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Package      FixPIMDB
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng (yuxing at uchicago.edu)

   Updated      Oct-01-2011
   Version      1.0
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include "fix_pimdb.h"
#include "universe.h"
#include "comm.h"
#include "force.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include <algorithm> 

using namespace LAMMPS_NS;
using namespace FixConst;

enum{PIMD,NMPIMD,CMD};

/* ---------------------------------------------------------------------- */

FixPIMDB::FixPIMDB(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  method     = PIMD;
  fmass      = 1.0;
  nhc_temp   = 298.15;
  nhc_nchain = 2;
  sp         = 1.0;
  nevery     = 100;
  E_kn = std::vector<double>((atom->nlocal * (atom->nlocal + 1) / 2),0.0);
  V = std::vector<double>((atom->nlocal + 1),0.0);

  for(int i=3; i<narg-1; i+=2)
  {
    if(strcmp(arg[i],"method")==0)
    {
      if(strcmp(arg[i+1],"pimd")==0) method=PIMD;
      else if(strcmp(arg[i+1],"nmpimd")==0) method=NMPIMD;
      else if(strcmp(arg[i+1],"cmd")==0) method=CMD;
      else error->universe_all(FLERR,"Unkown method parameter for fix pimd-B");
    }
    else if(strcmp(arg[i],"fmass")==0)
    {
      fmass = atof(arg[i+1]);
      if(fmass<0.0 || fmass>1.0) error->universe_all(FLERR,"Invalid fmass value for fix pimd-B");
    }
    else if(strcmp(arg[i],"sp")==0)
    {
      sp = atof(arg[i+1]);
      if(fmass<0.0) error->universe_all(FLERR,"Invalid sp value for fix pimd-B");
    }
    else if(strcmp(arg[i],"temp")==0)
    {
      nhc_temp = atof(arg[i+1]);
      if(nhc_temp<0.0) error->universe_all(FLERR,"Invalid temp value for fix pimd-B");
    }
    else if(strcmp(arg[i],"nhc")==0)
    {
      nhc_nchain = atoi(arg[i+1]);
      if(nhc_nchain<2) error->universe_all(FLERR,"Invalid nhc value for fix pimd-B");
    }
    else error->universe_all(arg[i],i+1,"Unkown keyword for fix pimd-B");
  }

  /* Initiation */

  max_nsend = 0;
  tag_send = NULL;
  buf_send = NULL;

  max_nlocal = 0;
  buf_recv = NULL;
  buf_beads = NULL;

  size_plan = 0;
  plan_send = plan_recv = NULL;

  M_x2xp = M_xp2x = M_f2fp = M_fp2f = NULL;
  lam = NULL;
  mode_index = NULL;

  mass = NULL;

  array_atom = NULL;
  nhc_eta = NULL;
  nhc_eta_dot = NULL;
  nhc_eta_dotdot = NULL;
  nhc_eta_mass = NULL;

  size_peratom_cols = 12 * nhc_nchain + 3;

  nhc_offset_one_1 = 3 * nhc_nchain;
  nhc_offset_one_2 = 3 * nhc_nchain +3;
  nhc_size_one_1 = sizeof(double) * nhc_offset_one_1;
  nhc_size_one_2 = sizeof(double) * nhc_offset_one_2;

  restart_peratom = 1;
  peratom_flag    = 1;
  peratom_freq    = 1;

  global_freq = 1;
  thermo_energy = 1;
  vector_flag = 1;
  size_vector = 2;
  extvector   = 1;
  comm_forward = 3;

  atom->add_callback(0); // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(1); // Call LAMMPS to re-assign restart-data for per-atom array

  grow_arrays(atom->nmax);

  // some initilizations

  nhc_ready = false;
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::init()
{
  if (atom->map_style == 0)
    error->all(FLERR,"Fix pimd-B requires an atom map, see atom_modify");

  if(universe->me==0 && screen) fprintf(screen,"Fix pimd-B initializing Path-Integral ...\n");

  // prepare the constants

  np = universe->nworlds;
  inverse_np = 1.0 / np;

  /* The first solution for the force constant, using SI units

  const double Boltzmann = 1.3806488E-23;    // SI unit: J/K
  const double Plank     = 6.6260755E-34;    // SI unit: m^2 kg / s

  double hbar = Plank / ( 2.0 * M_PI ) * sp;
  double beta = 1.0 / ( Boltzmann * input.nh_temp);

  // - P / ( beta^2 * hbar^2)   SI unit: s^-2
  double _fbond = -1.0 / (beta*beta*hbar*hbar) * input.nbeads;

  // convert the units: s^-2 -> (kcal/mol) / (g/mol) / (A^2)
  fbond = _fbond * 4.184E+26;

  */

  /* The current solution, using LAMMPS internal real units */

  const double Boltzmann = force->boltz;
  const double Plank     = force->hplanck;

  double hbar   = Plank / ( 2.0 * M_PI );
  double beta   = 1.0 / (Boltzmann * nhc_temp);
  double _fbond = 1.0 * np / (beta*beta*hbar*hbar) ;

  omega_np = sqrt(np) / (hbar * beta) * sqrt(force->mvv2e);
  fbond = - _fbond * force->mvv2e;

  if(universe->me==0)
    printf("Fix pimd-B -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  comm_init();

  mass = new double [atom->ntypes+1];

  if(method==CMD || method==NMPIMD) nmpimd_init();
  else for(int i=1; i<=atom->ntypes; i++) mass[i] = atom->mass[i] / np * fmass;

  if(!nhc_ready) nhc_init();
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::setup(int vflag)
{
  if(universe->me==0 && screen) fprintf(screen,"Setting up Path-Integral ...\n");

  post_force(vflag);
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::initial_integrate(int /*vflag*/)
{
  nhc_update_v();
  nhc_update_x();
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::final_integrate()
{
  nhc_update_v();
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::post_force(int /*flag*/)
{
  for(int i=0; i<atom->nlocal; i++) for(int j=0; j<3; j++) atom->f[i][j] /= np;

  comm_exec(atom->x);
  spring_force();

  if(method==CMD || method==NMPIMD)
  {
    /* forward comm for the force on ghost atoms */

    nmpimd_fill(atom->f);

    /* inter-partition comm */

    comm_exec(atom->f);

    /* normal-mode transform */

    nmpimd_transform(buf_beads, atom->f, M_f2fp[universe->iworld]);
  }
}

/* ----------------------------------------------------------------------
   Nose-Hoover Chains
------------------------------------------------------------------------- */

void FixPIMDB::nhc_init()
{
  double tau = 1.0 / omega_np;
  double KT  = force->boltz * nhc_temp;

  double mass0 = KT * tau * tau;
  int max = 3 * atom->nlocal;

  for(int i=0; i<max; i++)
  {
    for(int ichain=0; ichain<nhc_nchain; ichain++)
    {
      nhc_eta[i][ichain]        = 0.0;
      nhc_eta_dot[i][ichain]    = 0.0;
      nhc_eta_dot[i][ichain]    = 0.0;
      nhc_eta_dotdot[i][ichain] = 0.0;
      nhc_eta_mass[i][ichain]   = mass0;
      if((method==CMD || method==NMPIMD) && universe->iworld==0) ; else nhc_eta_mass[i][ichain]  *= fmass;
    }

    nhc_eta_dot[i][nhc_nchain]    = 0.0;

    for(int ichain=1; ichain<nhc_nchain; ichain++)
      nhc_eta_dotdot[i][ichain] = (nhc_eta_mass[i][ichain-1] * nhc_eta_dot[i][ichain-1]
        * nhc_eta_dot[i][ichain-1] * force->mvv2e - KT) / nhc_eta_mass[i][ichain];
  }

  // Zero NH acceleration for CMD

  if(method==CMD && universe->iworld==0) for(int i=0; i<max; i++)
    for(int ichain=0; ichain<nhc_nchain; ichain++) nhc_eta_dotdot[i][ichain] = 0.0;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::nhc_update_x()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;

  if(method==CMD || method==NMPIMD)
  {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);

    /* borrow the space of atom->f to store v in cartisian */

    v = atom->f;
    nmpimd_transform(buf_beads, v, M_xp2x[universe->iworld]);
  }

  for(int i=0; i<n; i++)
  {
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::nhc_update_v()
{
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for(int i=0; i<n; i++)
  {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

  t_sys = 0.0;
  if(method==CMD && universe->iworld==0) return;

  double expfac;
  int nmax = 3 * atom->nlocal;
  double KT = force->boltz * nhc_temp;
  double kecurrent, t_current;

  double dthalf = 0.5   * update->dt;
  double dt4    = 0.25  * update->dt;
  double dt8    = 0.125 * update->dt;

  for(int i=0; i<nmax; i++)
  {
    int iatm = i/3;
    int idim = i%3;

    double *vv = v[iatm];

    kecurrent = mass[type[iatm]] * vv[idim]* vv[idim] * force->mvv2e;
    t_current = kecurrent / force->boltz;

    double *eta = nhc_eta[i];
    double *eta_dot = nhc_eta_dot[i];
    double *eta_dotdot = nhc_eta_dotdot[i];

    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for(int ichain=nhc_nchain-1; ichain>0; ichain--)
    {
      expfac = exp(-dt8 * eta_dot[ichain+1]);
      eta_dot[ichain] *= expfac;
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    expfac = exp(-dt8 * eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    // Update particle velocities half-step

    double factor_eta = exp(-dthalf * eta_dot[0]);
    vv[idim] *= factor_eta;

    t_current *= (factor_eta * factor_eta);
    kecurrent = force->boltz * t_current;
    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for(int ichain=0; ichain<nhc_nchain; ichain++)
      eta[ichain] += dthalf * eta_dot[ichain];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    for(int ichain=1; ichain<nhc_nchain; ichain++)
    {
      expfac = exp(-dt8 * eta_dot[ichain+1]);
      eta_dot[ichain] *= expfac;
      eta_dotdot[ichain] = (nhc_eta_mass[i][ichain-1] * eta_dot[ichain-1] * eta_dot[ichain-1]
                           - KT) / nhc_eta_mass[i][ichain];
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    t_sys += t_current;
  }

  t_sys /= nmax;
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixPIMDB::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");
  memory->create(M_f2fp, np, np, "fix_feynman:M_f2fp");
  memory->create(M_fp2f, np, np, "fix_feynman:M_fp2f");

  lam = (double*) memory->smalloc(sizeof(double)*np, "FixPIMDB::lam");

  // Set up  eigenvalues

  lam[0] = 0.0;
  if(np%2==0) lam[np-1] = 4.0 * np;

  for(int i=2; i<=np/2; i++)
  {
    lam[2*i-3] = lam[2*i-2] = 2.0 * np * (1.0 - 1.0 *cos(2.0*M_PI*(i-1)/np));
  }

  // Set up eigenvectors for non-degenerated modes

  for(int i=0; i<np; i++)
  {
    M_x2xp[0][i] = 1.0 / np;
    if(np%2==0) M_x2xp[np-1][i] = 1.0 / np * pow(-1.0, i);
  }

  // Set up eigenvectors for degenerated modes

  for(int i=0; i<(np-1)/2; i++) for(int j=0; j<np; j++)
  {
    M_x2xp[2*i+1][j] =   sqrt(2.0) * cos ( 2.0 * M_PI * (i+1) * j / np) / np;
    M_x2xp[2*i+2][j] = - sqrt(2.0) * sin ( 2.0 * M_PI * (i+1) * j / np) / np;
  }

  // Set up Ut

  for(int i=0; i<np; i++)
    for(int j=0; j<np; j++)
    {
      M_xp2x[i][j] = M_x2xp[j][i] * np;
      M_f2fp[i][j] = M_x2xp[i][j] * np;
      M_fp2f[i][j] = M_xp2x[i][j];
    }

  // Set up masses

  int iworld = universe->iworld;

  for(int i=1; i<=atom->ntypes; i++)
  {
    mass[i] = atom->mass[i];

    if(iworld)
    {
      mass[i] *= lam[iworld];
      mass[i] *= fmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::nmpimd_fill(double **ptr)
{
  comm_ptr = ptr;
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::nmpimd_transform(double** src, double** des, double *vector)
{
  int n = atom->nlocal;
  int m = 0;

  for(int i=0; i<n; i++) for(int d=0; d<3; d++)
  {
    des[i][d] = 0.0;
    for(int j=0; j<np; j++) { des[i][d] += (src[j][m] * vector[j]); }
    m++;
  }
}

//dE_n^(k) is a function of k atoms (R_n-k+1,...,R_n) for a given n and k.
std::vector<double> FixPIMDB::Evaluate_dEkn_on_atom(const int n, const int k, const int atomnum)
{
  //dE_n^(k)(R_n-k+1,...,R_n) is a function of k atoms
  if (atomnum < n-k or atomnum > n-1 ) { return std::vector<double>(3, 0.0); }
  else {

    //bead is the bead number of current replica. bead = 0,...,np-1.
    int bead = universe->iworld;

    double **x = atom->x;
    double *_mass = atom->mass;
    int *type = atom->type;
    int nlocal = atom->nlocal;

    //xnext is a pointer to first element of buf_beads[x_next].
    //See in FixPIMDB::comm_init() for the definition of x_next.
    //x_next is basically (bead + 1) for bead in (0,...,np-2) and 0 for bead = np-1.
    //buf_beads[j] is a 1-D array of length 3*nlocal x0^j,y0^j,z0^j,...,x_(nlocal-1)^j,y_(nlocal-1)^j,z_(nlocal-1)^j.
    double *xnext = buf_beads[x_next];
    double *xlast = buf_beads[x_last];

    //omega^2, could use fbond instead?
    double omega_sq = omega_np * omega_np;

    //dE_n^(k)(R_n-k+1,...,R_n) is a function of k atoms
    //But derivative if for atom atomnum
    xnext += 3 * (atomnum);
    xlast += 3 * (atomnum);

    //np is total number of beads
    if (bead == np - 1 && k > 1){
      atomnum == n - 1 ? (xnext-= 3*(k - 1)) : (xnext += 3);
    }

    if (bead == 0 && k > 1){
      atomnum == n-k ? (xlast+= 3*(k - 1)) : (xlast -= 3);
    }

    //if (bead == np - 1 && k > 1) xnext += 3;
    //if (bead == 0 && k > 1) xlast -= 3;

    /*
    if(bead==3 && n==1 && k==1) {
      std::cout << "atom " << atomnum + 1 << ", bead" << bead + 1 << ": " << x[atomnum][0] << " " << x[atomnum][1] << " " << x[atomnum][2]
                << std::endl;
      std::cout << "next " << atomnum + 1 << ", bead" << bead + 1 << ": " << xnext[0] << " " << xnext[1] << " " << xnext[2]
                << std::endl;
    }*/

    std::vector<double> res(3);
/*
    std::cout<< "atom " << atomnum+1 << ", bead" << bead + 1 << ": " << x[atomnum][0] << " " << x[atomnum][1] << " " << x[atomnum][2] << std::endl;
    std::cout<< "next " << atomnum+1 << ", bead" << bead + 1 << ": " << xnext[0] << " " << xnext[1] << " " << xnext[2] << std::endl;
    std::cout<< "last " << atomnum+1 << ", bead" << bead + 1 << ": " << xlast[0] << " " << xlast[1] << " " << xlast[2] << std::endl;
*/
    double delx1 = xnext[0] - x[atomnum][0];
    double dely1 = xnext[1] - x[atomnum][1];
    double delz1 = xnext[2] - x[atomnum][2];
    domain->minimum_image(delx1, dely1, delz1);

    double delx2 = xlast[0] - x[atomnum][0];
    //std::cout<<  xnext[0] << std::endl;
    double dely2 = xlast[1] - x[atomnum][1];
    double delz2 = xlast[2] - x[atomnum][2];
    domain->minimum_image(delx2, dely2, delz2);

    double dx = -1.0*(delx1 + delx2);
    double dy = -1.0*(dely1 + dely2);
    double dz = -1.0*(delz1 + delz2);

    //std::cout << delx << " " <<dely << " " <<  delz << std::endl;
    //std::cout << _mass[type[i]] << " " << omega_sq << " " <<  delx*delx << std::endl;
    res.at(0) = _mass[type[atomnum]] * omega_sq * dx;
    res.at(1) = _mass[type[atomnum]] * omega_sq * dy;
    res.at(2) = _mass[type[atomnum]] * omega_sq * dz;

    //std::cout << bead << ": " << res.at(0) << " " << res.at(1) << " " << res.at(2) << std::endl;

    return res;
  }

}


/* ---------------------------------------------------------------------- */

//E_n^(k) is a function of k atoms (R_n-k+1,...,R_n) for a given n and k.
double FixPIMDB::Evaluate_Ekn(const int n, const int k)
{
  //bead is the bead number of current replica. bead = 0,...,np-1.
  int bead = universe->iworld;

  double **x = atom->x;
  double* _mass = atom->mass;
  int* type = atom->type;
  int nlocal = atom->nlocal;

  //xnext is a pointer to first element of buf_beads[x_next].
  //See in FixPIMDB::comm_init() for the definition of x_next.
  //x_next is basically (bead + 1) for bead in (0,...,np-2) and 0 for bead = np-1.
  //buf_beads[j] is a 1-D array of length 3*nlocal x0^j,y0^j,z0^j,...,x_(nlocal-1)^j,y_(nlocal-1)^j,z_(nlocal-1)^j.
  double* xnext = buf_beads[x_next];

  //omega^2, could use fbond instead?
  double omega_sq = omega_np*omega_np;

  //E_n^(k)(R_n-k+1,...,R_n) is a function of k atoms
  xnext += 3*(n-k);

  //np is total number of beads
  if(bead == np-1 && k > 1) xnext += 3;

  spring_energy = 0.0;
  for (int i = n-k; i < n ; ++i) {

    /*if(bead==3 && n==2 && k==2) {
        std::cout << "atom " << i + 1 << ", bead" << bead + 1 << ": " << x[i][0] << " " << x[i][1] << " " << x[i][2]
                  << std::endl;
        std::cout << "next " << i + 1 << ", bead" << bead + 1 << ": " << xnext[0] << " " << xnext[1] << " " << xnext[2]
                  << std::endl;
    }*/

    double delx = xnext[0] - x[i][0];
    //std::cout<<  xnext[0] << std::endl;
    double dely = xnext[1] - x[i][1];
    double delz = xnext[2] - x[i][2];

    domain->minimum_image(delx, dely, delz);

    if (bead == np - 1 && i == n - 2) {
      /*if(bead==3 && n==2 && k==2) {
          std::cout<<"I AM HERE"<<std::endl;
          std::cout << "next " << i + 1 << ", bead" << bead + 1 << ": " << xnext[0] << " " << xnext[1] << " " << xnext[2]
                    << std::endl;
      }*/
      xnext = buf_beads[x_next];

      /*if(bead==3 && n==2 && k==2) {
          std::cout<<"NOW I AM HERE"<<std::endl;
          std::cout << "next " << i + 1 << ", bead" << bead + 1 << ": " << xnext[0] << " " << xnext[1] << " " << xnext[2]
                    << std::endl;
      }*/
      //std::cout<<bead<<std::endl;
      //std::cout<<  xnext[0] << " " << xnext[1]<< " " << xnext[2] << std::endl;

      xnext += 3*(n - k);
    } else xnext += 3;

    //std::cout << delx << " " <<dely << " " <<  delz << std::endl;
    //std::cout << _mass[type[i]] << " " << omega_sq << " " <<  delx*delx << std::endl;
    spring_energy += 0.5*_mass[type[i]]*omega_sq*(delx*delx + dely*dely + delz*delz);

  }

  double energy_all = 0.0;
  double energy_local = spring_energy;
  //double energy_local = 0.0;
  //if(bead==0 && n==2 && k==2)
      //std::cout<< universe->iworld << " " << spring_energy <<" " << energy_all <<std::endl;

  //MPI_Allreduce(&spring_energy,&energy_local,1,MPI_DOUBLE,MPI_SUM,world);
  //MPI_Allreduce(&energy_local,&energy_all,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  MPI_Allreduce(MPI_IN_PLACE,&energy_local,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  if(isnan(spring_energy) || isnan(energy_all)){
    std::cout<< universe->iworld << " " << spring_energy <<" " << energy_all <<std::endl;
    exit(0);}

  //return energy_all;
  return energy_local;


}

std::vector<std::vector<double>>
FixPIMDB::Evaluate_dVBn(const std::vector<double> &V, const std::vector<double> &save_E_kn, const int n) {

  const double Boltzmann = force->boltz;
  double beta   = 1.0 / (Boltzmann * nhc_temp);
  int bead = universe->iworld;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  std::vector<std::vector<double>> dV_all(n, std::vector<double>(3,0.0));
  //std::cout<<dV_all.at(0).size()<<std::endl;

  for (int atomnum = 0; atomnum < nlocal; ++atomnum) {

      std::vector<std::vector<double>> dV(n+1, std::vector<double>(3,0.0));
      dV.at(0) = {0.0,0.0,0.0};

      for (int m = 1; m < n + 1; ++m) {

        std::vector<double> sig(3,0.0);

        if (atomnum > m-1) {
          dV.at(m) = {0.0,0.0,0.0};
        }else{

          //for (int k = 1; k < m + 1; ++k) {
          int count = m*(m-1)/2;

	  //!BH! I had to add 0.5 below in order to not get sigma which is zero or inf for large systems
          double Elongest = 0.5*(save_E_kn.at(m*(m-1)/2)+V.at(m-1));
	  //double Elongest = save_E_kn.at(m*(m-1)/2);
          //std::cout<<"Elongest: "<<Elongest<<std::endl;
/*
          for (int k = m; k > 0; --k) {

            std::vector<double> dE_kn(3,0.0);

            dE_kn = Evaluate_dEkn_on_atom(m,k,atomnum);

            sig.at(0) += (dE_kn.at(0) + dV.at(m - k).at(0)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));
            sig.at(1) += (dE_kn.at(1) + dV.at(m - k).at(1)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));
            sig.at(2) += (dE_kn.at(2) + dV.at(m - k).at(2)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));

            count++;


          }



          std::vector<double> num(3,0.0);
          num.at(0) = Elongest-1.0/beta*log(sig.at(0));
          num.at(1) = Elongest-1.0/beta*log(sig.at(1));
          num.at(2) = Elongest-1.0/beta*log(sig.at(2));

          //double  denom = -1.0/beta*log((double)m)+ V.at(m);


          //std::cout<<m<<" "<<beta<<" "<<V.at(m)<<std::endl;
          //std::cout<<"sig[0] is: " << sig[0] << " V.at(m) is: "<< V.at(m) <<std::endl;
          //dV.at(m).at(0) = exp(-beta*(num.at(0)-denom));
          //dV.at(m).at(1) = exp(-beta*(num.at(1)-denom));
          //dV.at(m).at(2) = exp(-beta*(num.at(2)-denom));
          double  sig_denom_m = (double)m*exp(-beta*V.at(m));
          dV.at(m).at(0) = sig.at(0) / sig_denom_m;
          dV.at(m).at(1) = sig.at(1) / sig_denom_m;
          dV.at(m).at(2) = sig.at(2) / sig_denom_m;
*/
          //std::cout<<"m: " << m<<" beta: "<<beta <<" num[0]: " <<num.at(0)<<" denom: "<<denom<<std::endl;
          //exit(1);

            for (int k = m; k > 0; --k) {
                std::vector<double> dE_kn(3,0.0);

                dE_kn = Evaluate_dEkn_on_atom(m,k,atomnum);
                /*
                if(bead==0 && atomnum==1) {
                    std::cout << "m: " << m <<  " k: " << k <<" Ekn:" << save_E_kn.at(count)* 2625.499638
                            <<" V_k-m:" << V.at(m - k)* 2625.499638 << std::endl;
                }*/

                sig.at(0) += (dE_kn.at(0) + dV.at(m - k).at(0)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));
                sig.at(1) += (dE_kn.at(1) + dV.at(m - k).at(1)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));
                sig.at(2) += (dE_kn.at(2) + dV.at(m - k).at(2)) * exp(-beta * (save_E_kn.at(count) + V.at(m - k)-Elongest));

                count++;

            }

            double  sig_denom_m = (double)m*exp(-beta*(V.at(m)-Elongest));
	    if(sig_denom_m ==0 || isnan(sig_denom_m) || isinf(sig_denom_m) || isnan(sig.at(0)) || isinf(sig.at(0))) {
	      if (universe->iworld ==0){
		std::cout << "m is: "<< m << " Elongest is: " << Elongest << " V.at(m-1) is " <<V.at(m-1)<< " beta is: " << beta << " sig_denom_m is: " <<sig_denom_m << std::endl;
	      }
	    }

            //std::cout<<m<<" "<<beta<<" "<<V.at(m)<<std::endl;
            //std::cout<<"sig[0] is: " <<sig.at(0)<<" sig_denom_m is: "<<sig_denom_m<<std::endl;
            dV.at(m).at(0) = sig.at(0) / sig_denom_m;
            dV.at(m).at(1) = sig.at(1) / sig_denom_m;
            dV.at(m).at(2) = sig.at(2) / sig_denom_m;

	    if(isinf(dV.at(m).at(0)) || isnan(dV.at(m).at(0))) {
	      if (universe->iworld ==0){
		std::cout << "sig_denom_m is: " << sig_denom_m << " Elongest is: " << Elongest
			  << " V.at(m) is " << V.at(m) << " beta is " << beta << std::endl;}
	      exit(0);
	    }


        }


      }

      /*if(bead==0 &&atomnum==0)
          std::cout <<"atom: " << atomnum+1 <<" bead: "<< bead+1 << " dV: " << dV.at(n).at(0)<<" "<< dV.at(n).at(1)<< " "<<dV.at(n).at(2) <<std::endl;
      */

      //std::cout<<"index: " <<(atomnum)*np + (bead)<<std::endl;
      dV_all.at((atomnum)).at(0) = dV.at(n).at(0);
      dV_all.at((atomnum)).at(1) = dV.at(n).at(1);
      dV_all.at((atomnum)).at(2) = dV.at(n).at(2);

      /*if(bead==0)
          std::cout <<"atom: " << atomnum+1 <<" bead: "<< bead+1 << " fbefore: " << f[atomnum][0]<<" "<< f[atomnum][1]<< " "<<f[atomnum][2] <<std::endl;
      */

      f[atomnum][0] -= dV.at(n).at(0);
      f[atomnum][1] -= dV.at(n).at(1);
      f[atomnum][2] -= dV.at(n).at(2);

      /*if(bead==0)
          std::cout <<"atom: " << atomnum+1 <<" bead: "<< bead+1 << " fafter: " << f[atomnum][0]<<" "<< f[atomnum][1]<< " "<<f[atomnum][2] <<std::endl;
        */
  }


  return dV_all;

}

std::vector<double> FixPIMDB::Evaluate_VBn(std::vector <double>& V, const int n)
{
  const double Boltzmann = force->boltz;
  double beta   = 1.0 / (Boltzmann * nhc_temp);
  std::vector<double> save_E_kn(n*(n+1)/2);

  int count = 0;
  for (int m = 1; m < n+1; ++m) {
    double sig_denom = 0.0;

    double Elongest=0.0;
    //for (int k = 1; k < m+1; ++k) {
    for (int k = m; k > 0; --k) {
          double E_kn;

          E_kn = Evaluate_Ekn(m,k);
          if(k==m){
            //!BH! I had to add 0.5 below in order to not get sigma which is zero or inf for large systems
	    Elongest = 0.5*(E_kn+V.at(m-1));
	    //Elongest = E_kn;
	    //Elongest = 0.5*E_kn;
	    //Elongest = 0.5*(std::max(E_kn,V.at(m-1)));
	  }
            //std::cout<<"Elongest: "<<Elongest<<std::endl;

          sig_denom += exp(-beta*(E_kn + V.at(m-k)-Elongest));
          //save_E_kn.push_back(E_kn);
          save_E_kn.at(count) = E_kn;
          count++;
	  if(sig_denom ==0 || isnan(sig_denom) || isinf(sig_denom)) {
	      if (universe->iworld ==0){
              std::cout << "m is: "<<m << " k is: " <<k << " E_kn is: " << E_kn << " V.at(m-k) is: " << V.at(m - k) << " Elongest is: " << Elongest
                        << " V.at(m-1) is " <<V.at(m-1)<< " beta is: " << beta << " sig_denom is: " <<sig_denom << std::endl;}
          }
    }

    V.at(m) = Elongest-1.0/beta*log(sig_denom / (double)m);
    if(isinf(V.at(m)) || isnan(V.at(m))) {
	if (universe->iworld ==0){
          std::cout << "sig_denom is: " << sig_denom << " Elongest is: " << Elongest
                    << std::endl;}
          exit(0);
    }
    //std::cout<< sig_denom << " " <<  log (sig_denom / (double)m) << " " << beta <<std::endl;
  }


  return save_E_kn;

}

/* ---------------------------------------------------------------------- */

void FixPIMDB::spring_force() {

    int nlocal = atom->nlocal;

    V.at(0) = 0.0;
    std::vector<std::vector<double>> dV(nlocal * universe->nworlds, std::vector<double>(3, 0.0));

    //if(universe->iworld ==0)
    //std::cout<<"replica\tn\tt_V[microsec]\tt_dV[microsec]"<<std::endl;

    //auto start = std::chrono::high_resolution_clock::now();
    E_kn = Evaluate_VBn(V, nlocal);

    //auto elapsed = std::chrono::high_resolution_clock::now() - start;
    //long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    ////std::cout<<nlocal<<"\t"<<microseconds<<";"<<std::endl;

    //auto start_d = std::chrono::high_resolution_clock::now();
    dV = Evaluate_dVBn(V,E_kn,nlocal);
    //auto elapsed_d = std::chrono::high_resolution_clock::now() - start_d;
    //long long microseconds_d = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_d).count();
    ////std::cout<<nlocal<<"\t"<<microseconds_d<<";"<<std::endl;

    //std::cout<<universe->iworld <<"\t"<<nlocal<<"\t"<<microseconds<<"\t"<<microseconds_d<<";"<<std::endl;


    //PRINTING FOR DEBUGGING/
/*
    if(universe->iworld ==0) {
      //std::cout<<nlocal<<"\t"<<microseconds<<"\t"<<microseconds_d<<";"<<std::endl;
      for (double val : V)
        std::cout << val << std::endl;
      std::cout << "" << std::endl;
      for (double val : E_kn)
        std::cout << val << std::endl;
      std::cout << "" << std::endl;

      for(std::vector<double> vec: dV) {
          std::cout << vec[0] <<" " <<vec[1] <<" " <<vec[2] <<std::endl;
          std::cout << std::endl;
      }

    }
*/

    /*
    int atomnum =1;
    int beadnum =3;
    if (universe->iworld ==beadnum) {
      for (int m = 1; m < atom->nlocal + 1; ++m) {
        for (int k = 1; k < m + 1; ++k) {
          std::vector<double> dE = Evaluate_dEkn_on_atom(m, k, atomnum);
          std::cout <<"atom " << atomnum+1 << " bead " << universe->iworld +1 << " n " << m << " k " << k << " " << dE.at(0) << " " << dE.at(1) << " " << dE.at(2) << std::endl;
        }
      }
    }*/

    //double E = Evaluate_Ekn(3,3);
    //std::vector<double> dE = Evaluate_dEkn_on_atom(3,3,0);

}

/* ----------------------------------------------------------------------
 OLD spring_force!!!!!!!!!!!!!!!!!!!!!!
 ---------------------------------------------------------------------- */
/*
void FixPIMDB::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double* _mass = atom->mass;
  int* type = atom->type;
  int nlocal = atom->nlocal;

  double* xlast = buf_beads[x_last];
  double* xnext = buf_beads[x_next];

  for(int i=0; i<nlocal; i++)
  {
    double delx1 = xlast[0] - x[i][0];
    double dely1 = xlast[1] - x[i][1];
    double delz1 = xlast[2] - x[i][2];
    xlast += 3;
    domain->minimum_image(delx1, dely1, delz1);

    double delx2 = xnext[0] - x[i][0];
    double dely2 = xnext[1] - x[i][1];
    double delz2 = xnext[2] - x[i][2];
    xnext += 3;
    domain->minimum_image(delx2, dely2, delz2);

    double ff = fbond * _mass[type[i]];

    double dx = delx1+delx2;
    double dy = dely1+dely2;
    double dz = delz1+delz2;

    f[i][0] -= (dx) * ff;
    f[i][1] -= (dy) * ff;
    f[i][2] -= (dz) * ff;

    spring_energy += (dx*dx+dy*dy+dz*dz);
  }
}
*/

//FOR PRINTING ENERGIES AND POTENTIALS FOR PIMD-B
void FixPIMDB::end_of_step() {

    if (universe->iworld == 0) {
    //std::cout << "E1\tE12\tE2\tE123\tE23\tE3\tVB0\tVB1\tVB2\tVB3" <<std::endl;
    //#! FIELDS time E1 E2 E12 E3 E23 E123 VB1 VB2 VB3 E_ox3 Vr.bias
      std::ofstream myfile;
      myfile.open ("pimdb.log", std::ios::out | std::ios::app);

      for (double val: E_kn)
        myfile << val * 2625.499638 << " ";
      for (double val: V)
        myfile << val * 2625.499638 << " ";
      myfile << std::endl;

      myfile.close();

    }

}
/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixPIMDB::comm_init()
{
  if(size_plan)
  {
    delete [] plan_send;
    delete [] plan_recv;
  }

  if(method == PIMD)
  {
    size_plan = 2;
    plan_send = new int [2];
    plan_recv = new int [2];
    mode_index = new int [2];

    int rank_last = universe->me - comm->nprocs;
    int rank_next = universe->me + comm->nprocs;
    if(rank_last<0) rank_last += universe->nprocs;
    if(rank_next>=universe->nprocs) rank_next -= universe->nprocs;

    plan_send[0] = rank_next; plan_send[1] = rank_last;
    plan_recv[0] = rank_last; plan_recv[1] = rank_next;

    mode_index[0] = 0; mode_index[1] = 1;
    x_last = 1; x_next = 0;
  }
  else
  {
    size_plan = np - 1;
    plan_send = new int [size_plan];
    plan_recv = new int [size_plan];
    mode_index = new int [size_plan];

    for(int i=0; i<size_plan; i++)
    {
      plan_send[i] = universe->me + comm->nprocs * (i+1);
      if(plan_send[i]>=universe->nprocs) plan_send[i] -= universe->nprocs;

      plan_recv[i] = universe->me - comm->nprocs * (i+1);
      if(plan_recv[i]<0) plan_recv[i] += universe->nprocs;

      mode_index[i]=(universe->iworld+i+1)%(universe->nworlds);
    }

    x_next = (universe->iworld+1+universe->nworlds)%(universe->nworlds);
    x_last = (universe->iworld-1+universe->nworlds)%(universe->nworlds);
  }

  if(buf_beads)
  {
    for(int i=0; i<np; i++) if(buf_beads[i]) delete [] buf_beads[i];
    delete [] buf_beads;
  }

  buf_beads = new double* [np];
  for(int i=0; i<np; i++) buf_beads[i] = NULL;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::comm_exec(double **ptr)
{
  int nlocal = atom->nlocal;

  if(nlocal > max_nlocal)
  {
    max_nlocal = nlocal+200;
    int size = sizeof(double) * max_nlocal * 3;
    buf_recv = (double*) memory->srealloc(buf_recv, size, "FixPIMDB:x_recv");

    for(int i=0; i<np; i++)
      buf_beads[i] = (double*) memory->srealloc(buf_beads[i], size, "FixPIMDB:x_beads[i]");
  }

  // copy local positions

  memcpy(buf_beads[universe->iworld], &(ptr[0][0]), sizeof(double)*nlocal*3);

  // go over comm plans

  for(int iplan = 0; iplan<size_plan; iplan++)
  {
    // sendrecv nlocal

    int nsend;

    MPI_Sendrecv( &(nlocal), 1, MPI_INT, plan_send[iplan], 0,
                  &(nsend),  1, MPI_INT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // allocate arrays

    if(nsend > max_nsend)
    {
      max_nsend = nsend+200;
      tag_send = (tagint*) memory->srealloc(tag_send, sizeof(tagint)*max_nsend, "FixPIMDB:tag_send");
      buf_send = (double*) memory->srealloc(buf_send, sizeof(double)*max_nsend*3, "FixPIMDB:x_send");
    }

    // send tags

    MPI_Sendrecv( atom->tag, nlocal, MPI_LMP_TAGINT, plan_send[iplan], 0,
                  tag_send,  nsend,  MPI_LMP_TAGINT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions

    double *wrap_ptr = buf_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tag_send[i]);

      if(index<0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tag_send[i],
                universe->iworld, comm->me, plan_recv[iplan],
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }

      memcpy(wrap_ptr, ptr[index], ncpy);
      wrap_ptr += 3;
    }

    // sendrecv x

    MPI_Sendrecv( buf_send, nsend*3,  MPI_DOUBLE, plan_recv[iplan], 0,
                  buf_recv, nlocal*3, MPI_DOUBLE, plan_send[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // copy x

    memcpy(buf_beads[mode_index[iplan]], buf_recv, sizeof(double)*nlocal*3);
  }
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::pack_forward_comm(int n, int *list, double *buf,
                             int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = comm_ptr[j][0];
    buf[m++] = comm_ptr[j][1];
    buf[m++] = comm_ptr[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    comm_ptr[i][0] = buf[m++];
    comm_ptr[i][1] = buf[m++];
    comm_ptr[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Memory operations
------------------------------------------------------------------------- */

double FixPIMDB::memory_usage()
{
  double bytes = 0;
  bytes = atom->nmax * size_peratom_cols * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::grow_arrays(int nmax)
{
  if (nmax==0) return;
  int count = nmax*3;

  memory->grow(array_atom, nmax, size_peratom_cols, "FixPIMDB::array_atom");
  memory->grow(nhc_eta,        count, nhc_nchain,   "FixPIMDB::nh_eta");
  memory->grow(nhc_eta_dot,    count, nhc_nchain+1, "FixPIMDB::nh_eta_dot");
  memory->grow(nhc_eta_dotdot, count, nhc_nchain,   "FixPIMDB::nh_eta_dotdot");
  memory->grow(nhc_eta_mass,   count, nhc_nchain,   "FixPIMDB::nh_eta_mass");
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::copy_arrays(int i, int j, int /*delflag*/)
{
  int i_pos = i*3;
  int j_pos = j*3;

  memcpy(nhc_eta       [j_pos], nhc_eta       [i_pos], nhc_size_one_1);
  memcpy(nhc_eta_dot   [j_pos], nhc_eta_dot   [i_pos], nhc_size_one_2);
  memcpy(nhc_eta_dotdot[j_pos], nhc_eta_dotdot[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_mass  [j_pos], nhc_eta_mass  [i_pos], nhc_size_one_1);
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::pack_exchange(int i, double *buf)
{
  int offset=0;
  int pos = i * 3;

  memcpy(buf+offset, nhc_eta[pos],        nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_dot[pos],    nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(buf+offset, nhc_eta_dotdot[pos], nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_mass[pos],   nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::unpack_exchange(int nlocal, double *buf)
{
  int offset=0;
  int pos = nlocal*3;

  memcpy(nhc_eta[pos],        buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos],    buf+offset, nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos],   buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::pack_restart(int i, double *buf)
{
  int offset=0;
  int pos = i * 3;
  buf[offset++] = size_peratom_cols+1;

  memcpy(buf+offset, nhc_eta[pos],        nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_dot[pos],    nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(buf+offset, nhc_eta_dotdot[pos], nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_mass[pos],   nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

void FixPIMDB::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i=0; i<nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  int pos = nlocal * 3;

  memcpy(nhc_eta[pos],        extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos],    extra[nlocal]+m, nhc_size_one_2); m += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos],   extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::maxsize_restart()
{
  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

int FixPIMDB::size_restart(int /*nlocal*/)
{
  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

double FixPIMDB::compute_vector(int n)
{
  if(n==0) { return spring_energy; }
  if(n==1) { return t_sys; }
  return 0.0;
}


