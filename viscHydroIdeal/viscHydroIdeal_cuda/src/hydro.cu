#include "hydro.h"

#define del_eps 1e-7

hydro::hydro(EoS* _eos, grid* _f,idb* _IDB,double _tau0, double _dt, cnvrt* _CN)
{
 eos = _eos;
 f = _f;
 CN = _CN;
 tau = _tau0;
 IDB = _IDB ;
 dt =_dt;
}

hydro::~hydro()
{
}


void hydro::set_dtau(double deltaTau)
{
  dt = deltaTau;
  if (dt > IDB->dx / 2. ||
      dt > IDB->dy / 2. /*|| dt > tau*IDB->deta */) {
    cout << "too big delta_tau " << dt << "  " << IDB->dx << "  " << IDB->dy
	 << "  " << tau * IDB->deta << endl;
    exit(1);
  }
}


void hydro::evolve()
{

  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx(); ix++)
	{      
	  cell *c = f->get_cell(ix, iy, iz );
	  c->save_Q_prev();
	  c->clear_flux();
	}

  // X-direction flux
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx()-1; ix++)
	{
	  hlle_flux(f->get_cell(ix,iy, iz), f->get_cell(ix+1,iy,iz), X_, PREDICT);
	  // cout<<iy<<"\t"<<iz<<"\t"<<ix<<endl;
	}
  
  // Y-direction flux
  for(int iz = 0; iz<f->get_neta(); iz++)
    for(int ix = 0; ix<f->get_nx(); ix++)
      for(int iy = 0; iy<f->get_ny()-1; iy++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy+1,iz), Y_, PREDICT);
	}
  
if(f->get_neta() > 1) // don,t calculate z-flux in 2+1D hydro
{  
  // Z-direction flux
  for(int ix = 0; ix<f->get_nx(); ix++)
    for(int iy = 0; iy<f->get_ny(); iy++)
      for(int iz = 0; iz<f->get_neta() -1; iz++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy,iz+1), Z_, PREDICT);
	}
}  
  
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix< f->get_nx(); ix++)
        {
	  cell *c = f->get_cell(ix, iy,iz);
	  sourcestep( PREDICT,  ix,  iy, iz,  tau);
	  c->update_Q_to_Qh_by_flux();
	  c->clear_flux();
	}
  
  // X-direction flux
  for(int iy = 0; iy<f->get_ny(); iy++)
    for(int iz = 0; iz<f->get_neta(); iz++)
      for(int ix = 0; ix<f->get_nx()-1; ix++)
	{
	  hlle_flux(f->get_cell(ix,iy, iz), f->get_cell(ix+1,iy,iz), X_, CORRECT);
	}
  
  
  // Y-direction flux
  for(int iz = 0; iz<f->get_neta(); iz++)
    for(int ix = 0; ix<f->get_nx(); ix++)
      for(int iy = 0; iy<f->get_ny()-1; iy++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy+1,iz), Y_, CORRECT);
	}
  
if(f->get_neta()> 1) // don,t calculate z-flux in 2+1D hydro
{  
  // Z-direction flux
  for(int ix = 0; ix<f->get_nx(); ix++)
    for(int iy = 0; iy<f->get_ny(); iy++)
      for(int iz = 0; iz<f->get_neta()-1; iz++)
	{
	  hlle_flux(f->get_cell(ix,iy,iz), f->get_cell(ix,iy,iz+1), Z_, CORRECT);
	}
  
}  
  for(int iy = 0; iy <f->get_ny(); iy++)
    for(int iz = 0; iz< f->get_neta(); iz++)
      for(int ix = 0; ix< f->get_nx(); ix++)
	{
	  cell *c = f->get_cell(ix, iy,iz); 
	  sourcestep( CORRECT,  ix,  iy, iz, tau);
	  c->update_by_flux();
	  c->clear_flux();	  
	}
  
  tau = tau+dt;                 // tau increased
  f->correct_imaginary_cells(); //boundary condition
  
   
}




void hydro::hlle_flux(cell* left, cell* right, int direction, int mode)
{
  const double dta = mode == 0 ? dt/2. : dt;
  double el,vxl,vyl,vzl,nbl,nql,nsl,pl;
  double er,vxr,vyr,vzr,nbr,nqr,nsr,pr;
  double El,Er;
  double Utl=0.,Uxl=0.,Uyl=0.,Uzl=0.,Ubl=0.,Uql=0.,Usl=0.;
  double Utr=0.,Uxr=0.,Uyr=0.,Uzr=0.,Ubr=0.,Uqr=0.,Usr=0.;
  double Ftl=0.,Fxl=0.,Fyl=0.,Fzl =0,Fbl=0., Fql =0, Fsl =0;
  double Ftr=0.,Fxr=0.,Fyr=0.,Fzr =0,Fbr=0., Fqr =0, Fsr =0;
  double csb,vb,bl=0.,br=0.;
  double flux[7]={0.0};
  double tauFactor;  // fluxes are also multiplied by tau
  double dx =0;
  
  if(mode == PREDICT)
    {
      left->get_right_var(eos, tau, el, pl, nbl,nql,nsl, vxl, vyl, vzl, direction);
      right->get_left_var(eos, tau, er, pr, nbr,nqr,nsr, vxr, vyr, vzr, direction);
      El = (el+pl)*(1.0/(1.0-vxl*vxl-vyl*vyl-vzl*vzl));
      Er = (er+pr)*(1.0/(1.0-vxr*vxr-vyr*vyr-vzr*vzr));
      tauFactor = tau ;
    }
  else
    {
      left->get_right_varH(eos, tau, el, pl, nbl,nql,nsl, vxl, vyl, vzl, direction );
      right->get_left_varH( eos, tau, er, pr, nbr,nqr,nsr, vxr, vyr, vzr, direction);
      El = (el+pl)*(1.0/(1.0-vxl*vxl-vyl*vyl-vzl*vzl));
      Er = (er+pr)*(1.0/(1.0-vxr*vxr-vyr*vyr-vzr*vzr));
      tauFactor = tau +(0.5*dt);
    }
  
  if (el < 0.)
    {
      el = 0.;
      pl = 0.;
    }
  if (er < 0.)
    {
      er = 0.;
      pr = 0.;
    }
  
  if (el < del_eps && er < del_eps ) return;  // *1) no flux calculation if both sides are empty cells
  
  double gammal = 1.0 / sqrt(1 - vxl * vxl - vyl * vyl - vzl*vzl );
  double gammar = 1.0 / sqrt(1 - vxr * vxr - vyr * vyr - vzr*vzr );
  Utl = gammal * gammal * (el + pl) - pl;
  Uxl = gammal * gammal * (el + pl) * vxl;
  Uyl = gammal * gammal * (el + pl) * vyl;
  Uzl = gammal * gammal * (el + pl) * vzl;
  Ubl = gammal * nbl;
  Uql = gammal * nql;
  Usl = gammal * nsl;
  
  
  Utr = gammar * gammar * (er + pr) - pr;
  Uxr = gammar * gammar * (er + pr) * vxr;
  Uyr = gammar * gammar * (er + pr) * vyr;
  Uzr = gammar * gammar * (er + pr) * vzr;
  Ubr = gammar * nbr;
  Uqr = gammar * nqr;
  Usr = gammar * nsr;
  
  
  if(direction == X_)
    {
      Ftl = Utl * vxl + pl * vxl;
      Fxl = Uxl * vxl + pl;
      Fyl = Uyl * vxl;
      Fzl = Uzl * vxl;
      Fbl = Ubl * vxl;
      Fql = Uql * vxl;
      Fsl = Usl * vxl;
      
      Ftr = Utr * vxr + pr * vxr;
      Fxr = Uxr * vxr + pr;
      Fyr = Uyr * vxr;
      Fzr = Uzr * vxr;
      Fbr = Ubr * vxr;
      Fqr = Uqr * vxr;
      Fsr = Usr * vxr;
      
      
      // for the case of constant c_s only
      csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vxl - vxr, 2));

      vb = (sqrt(El) * vxl + sqrt(Er) * vxr) / (sqrt(El) + sqrt(Er));
      bl = min(0., min((vb - csb) / (1 - vb * csb),
		       (vxl - eos->cs_()) / (1 - vxl * eos->cs_())));
      br = max(0., max((vb + csb) / (1 + vb * csb),
		       (vxr + eos->cs_()) / (1 + vxr * eos->cs_())));
      
      dx = f->get_dx();
      
      
      if (el == 0.) bl = -1.; 
      if (er == 0.) br = 1.;
    }
  
  
  
  if(direction == Y_)
    {
      Ftl = Utl * vyl + pl * vyl;
      Fxl = Uxl * vyl ;
      Fyl = Uyl * vyl + pl;
      Fzl = Uzl * vyl;
      Fbl = Ubl * vyl;
      Fql = Uql * vyl;
      Fsl = Usl * vyl;
      
      Ftr = Utr * vyr + pr * vyr;
      Fxr = Uxr * vyr ;
      Fyr = Uyr * vyr + pr;
      Fzr = Uzr * vyr;
      Fbr = Ubr * vyr;
      Fqr = Uqr * vyr;
      Fsr = Usr * vyr;
      
      


      // for the case of constant c_s only
      
     csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vyl - vyr, 2));
     vb = (sqrt(El) * vyl + sqrt(Er) * vyr) / (sqrt(El) + sqrt(Er));
     bl = min(0., min((vb - csb) / (1 - vb * csb),
		      (vyl - eos->cs_()) / (1 - vyl * eos->cs_())));
     br = max(0., max((vb + csb) / (1 + vb * csb),
                   (vyr + eos->cs_()) / (1 + vyr * eos->cs_())));

   
     dx = f->get_dy();
     
     
     if (el == 0.) bl = -1.;
     if (er == 0.) br = 1.;
    }
  
  
  if(direction == Z_)
    {
      double tau1 =tauFactor;
      Ftl = Utl * vzl/tau1 + pl * vzl/tau1;
      Fxl = Uxl * vzl/tau1 ;
      Fyl = Uyl * vzl/tau1;
      Fzl = Uzl * vzl/tau1 + pl/tau1;
      Fbl = Ubl * vzl/tau1;
      Fql = Uql * vzl/tau1;
      Fsl = Usl * vzl/tau1;
      
      Ftr = Utr * vzr/tau1 + pr * vzr/tau1;
      Fxr = Uxr * vzr/tau1 ;
      Fyr = Uyr * vzr/tau1 ;
      Fzr = Uzr * vzr/tau1 + pr/tau1 ;
      Fbr = Ubr * vzr/tau1;
      Fqr = Uqr * vzr/tau1;
      Fsr = Usr * vzr/tau1;
      
      
      
    // for the case of constant c_s only
     
     csb = sqrt(eos->cs2_() +
		 0.5 * sqrt(El * Er) / pow(sqrt(El) + sqrt(Er), 2) *
                 pow(vzl - vzr, 2));

     vb = (sqrt(El) * vzl + sqrt(Er) * vzr) / (sqrt(El) + sqrt(Er));
     bl = (1.0/tau)*min(0., min((vb - csb) / (1 - vb * csb),
				(vzl - eos->cs_()) / (1 - vzl * eos->cs_())));
     br = (1.0/tau)*max(0., max((vb + csb) / (1 + vb * csb),
				 (vzr + eos->cs_()) / (1 + vzr * eos->cs_())));

    


     dx =f->get_deta();
     
     
    if (el == 0.) bl = -1./tau; 
    if (er == 0.) br = 1./tau;
    }
  

  if(bl == 0. && br == 0.) return;
  
  flux[T_] = tauFactor * dta / dx *
    (-bl * br * (Utl - Utr) + br * Ftl - bl * Ftr) / (-bl + br);
  flux[X_] = tauFactor * dta / dx *
    (-bl * br * (Uxl - Uxr) + br * Fxl - bl * Fxr) / (-bl + br);
  flux[Y_] = tauFactor * dta / dx *
    (-bl * br * (Uyl - Uyr) + br * Fyl - bl * Fyr) / (-bl + br);
  flux[Z_] = tauFactor * dta / dx *
    (-bl * br * (Uzl - Uzr) + br * Fzl - bl * Fzr) / (-bl + br);
  flux[NB_] = tauFactor * dta / dx *
    (-bl * br * (Ubl - Ubr) + br * Fbl - bl * Fbr) / (-bl + br);
  flux[NQ_] = tauFactor * dta / dx *
    (-bl * br * (Uql - Uqr) + br * Fql - bl * Fqr) / (-bl + br);
  flux[NS_] = tauFactor * dta / dx *
    (-bl * br * (Usl - Usr) + br * Fsl - bl * Fsr) / (-bl + br);
  
  
  left->add_flux(-flux[T_], -flux[X_], -flux[Y_], -flux[Z_], -flux[NB_],
		-flux[NQ_], -flux[NS_]);
  right->add_flux(flux[T_], flux[X_], flux[Y_], flux[Z_], flux[NB_], flux[NQ_],
		 flux[NS_]);

}




void hydro::sourcestep(int mode, int ix, int iy, int iz, double _tau)
{ 
  double S[7];
  double e,p,vx,vy,vz,nb,nq,ns;
  if( mode == PREDICT )
    {
      double _dt  = 0.5*dt; // (1/2)*dtau -> for predictor step 
      f->get_cell(ix,iy,iz)->get_Q(S);
      for(int i =0; i<7; i++){S[i] = S[i]/_tau; }
      CN->CALC_2_LRF(eos,  S, e, p, nb, nq, ns, vx, vy, vz);
      if (e < del_eps ) return;
      f->get_cell(ix,iy,iz)->add_flux( (-S[T_] * vz * vz - p * (1. + vz * vz))*_dt, 0.0, 0.0, -S[Z_]*_dt, 0.0,0.0, 0.0);  
    }
  else
    {
      _tau = _tau+0.5*dt;
      double _dt  = dt;
      f->get_cell(ix,iy,iz)->get_Qh(S);
      for(int i =0; i<7; i++){S[i] = S[i]/_tau; }
      CN->CALC_2_LRF(eos,  S, e, p, nb, nq, ns, vx, vy, vz);
      if (e < del_eps) return;
      f->get_cell(ix,iy,iz)->add_flux(  (-S[T_] * vz * vz - p * (1. + vz * vz))*_dt, 0.0, 0.0, (-S[Z_])*_dt, 0.0,0.0, 0.0);     
    }
  
}








