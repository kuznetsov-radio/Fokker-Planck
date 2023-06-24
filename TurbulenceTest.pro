pro TurbulenceTest
 restore, 'Fluxtube II.sav' ;restore the GX Simulator loop model
 L=6.345d7 ;loop length [m] - not saved in the .sav file
 T_loop=2.5d7 ;plasma temperature in the loop [K] - not saved in the .sav file
 
 smin=min(s)-(1d0-max(s)+min(s))/2 ;coordinate of the 1st (left) footpoint [relative units]
 
 s*=L ;convert the distance to absolute units [m]
 smin*=L ;convert the coordinate of the 1st footpoint to absolute units [m]
 n_th*=1d6 ;convert the thermal plasma density to SI units
 
 restore, 'injection.sav' ;restore the injection model
 k=5 ;location of the injection peak
  
 ;parameters of the injection function at the injection peak:
 nb=nb[k]*1d6 ;convert to SI units
 q_0=q_0[k]
 q_2=q_2[k]

 libname='FokkerPlankSolver.dll' ;Stupin's library name
 
 ;fundamental constants (SI units):
 qe=1.602176634d-19
 me=9.1093837015d-31
 c=299792458d0
 eps0=8.8541878128d-12
 kB=1.380649d-23
 
 ;model dimensions:
 n_E=60L
 n_mu=40L ;must be even
 n_z=1000L ;note that the actual number of nodes will be 2*n_z+1
 nz_loop_param=n_elements(s) ;number of voxels along the loop in the original GX Simulator model (restored)
 
 switch_turb=0L ;turbulence key
 switch_reverse_current=0L ;reverse current key
 general_param=[n_E, n_mu, n_z, nz_loop_param, switch_reverse_current, switch_turb] ;integer parameters of the model
 
 Emin=E0 ;minimum energy [MeV] (restored from the injection model)
 Emax=E1 ;maximum energy [MeV] (restored from the injection model)
  
 zmax=L ;maximum z coordinate = length of the loop [m]
 r_gridz=0.25*zmax ;parameter of the non-uniform z grid - 0.25*zmax is sufficient
                   ;means that exp(-zmax/2/r_gridz)=0.135 of all grid points are located outside the loop
  
 tmax=5d0      ;end time of the calculations [s]
 dt0=0.01      ;initial integration time step (at t=0) [s]
 dt_max=0.01   ;maximum allowed integration time step [s]
 dt_record=0.1 ;time interval between consequent saved states [s]
 
 gamma=0.1d0 ;parameter that controls the integration accuracy and integration step: 
             ;relative change of the distrubution function at one step cannot exceed gamma
  
 ;parameters of the non-uniform energy grid (chosen to provide the energy range from Emin to Emax - do not change):
 k_minE=1d0-exp(alog(1d0*n_E)/(1d0-Emax/Emin))
 r_gridE=-Emin/alog(1d0-k_minE)
 
 grid_param=[Emin, zmax, tmax, dt0, dt_max, gamma, r_gridE, r_gridz] ;floating-point parameters of the model
 
 n_inf=1d20 ;chromospheric (outside the loop) plasma density [m^{-3}]
 T_inf=10d3 ;chromospheric plasma temperature [K]
 log_inf=(T_inf lt 1.4e5) ? 9.1-0.5*alog(n_inf/1d6)+1.5*alog(T_inf) : $
                            15.1-0.5*alog(n_inf/1d6)+alog(T_inf) ;Coulomb logarithm
 conductivity_inf=3d0/4/sqrt(2d0*!dpi)*(4d0*!dpi*eps0)^2*(kB*T_inf)^1.5/qe^2/sqrt(me)/log_inf ;conductivity [SI units]
 inf_param=[n_inf, log_inf, conductivity_inf] ;plasma parameters in the chromosphere
 
 z_loop=s-smin ;input z grid (shifted to start from zero) [m]
 B_loop=double(B) ;input magnetic field along the loop (units do not matter)
 n_loop=double(n_th) ;input plasma density along the loop [m^{-3}]
 coulomb_log=dblarr(n_elements(s))
 for i=0, n_elements(s)-1 do coulomb_log[i]=(T_loop lt 1.4e5) ? 9.1-0.5*alog(n_loop[i]/1d6)+1.5*alog(T_loop) : $
                                                                15.1-0.5*alog(n_loop[i]/1d6)+alog(T_loop) ;Coulomb logarithm
 conductivity=3d0/4/sqrt(2d0*!dpi)*(4d0*!dpi*eps0)^2*(kB*T_loop)^1.5/qe^2/sqrt(me)/coulomb_log ;conductivity [SI units]
 loop_param=[transpose(z_loop), $
             transpose(B_loop), $
             transpose(n_loop), $
             transpose(coulomb_log), $
             transpose(conductivity)] ;loop parameters
 loop_param=transpose(loop_param) ;transpose the array for consistency with the DLL definitions

 ;turbulence parameters, lambda=lambda_0*(E/E_0)^{-a}:
 a_turb=1d0 ;power-law index; Musset et al. (2018): a_turb=0.95
 lambda_turb0=3d6 ;mean free path at E=Emin [m]; Musset et al. (2018): lambda_turb0=3.34d6 m
 E_turb=-r_gridE*alog(1d0-k_minE-(1d0-k_minE)*dindgen(n_E+1)/n_E) ;energy grid - same as in the DLL (do not change)
 turb_param=lambda_turb0*(E_turb/Emin)^(-a_turb) ;resulting mean free path

 ;initialize the DLL using the above-mentioned parameters: 
 res=call_external(libname, 'CREATE_SOLVER', general_param, grid_param, inf_param, loop_param, turb_param)
 
 ;reserve space for the E, mu, z, etc. grids, create the irregular grids, and extrapolate the input data to these grids:
 E_grid=dblarr(n_E+1)
 mu_grid=dblarr(n_mu+1)
 z_grid=dblarr(2*n_z+1)
 dlnB_out=dblarr(2*n_z+1)
 n_loop_out=dblarr(2*n_z+1)
 res=call_external(libname, 'GET_GRIDS', E_grid, mu_grid, z_grid, dlnB_out, n_loop_out)

 s_grid=(z_grid+smin)/L ;convert the z grid to relative units - this is needed to compute the distribution function
 
 ;create the initial electron distribution function:
 f=dblarr(n_E+1, n_mu+1, 2*n_z+1) ;reserve space
 Fx=(gam-1)*(E_grid/E0)^(-gam)/E0 ;energy distribution (the restored parameters gam and E0 are used)
 Fmu=0.5 ;pitch-angle distribution
 Fz=E0*exp(-(q_0*(s_grid+q_2))^2-(q_0*(s_grid+q_2))^4) ;spatial distribution
 dt_inj=1d0 ;duration of injection (injection rate is assumed constant and equal Fx*Fmu*Fz)
 for ii=0, n_E do for jj=0, n_mu do for kk=0, 2*n_z do f[ii, jj, kk]=nb*Fx[ii]*Fmu*Fz[kk]*dt_inj ;locally injected distribution
; for ii=0, n_E do for jj=0, n_mu do f[ii, jj, *]+=mean(f[ii, jj, *]) ;add a uniform background
 f=transpose(f) ;transpose the array for consistency with the DLL definitions          

 ;other input/output parameters: 
 J=dblarr(2*n_z+1) ;return current (input/output), is assumed to be zero
 n_fast=dblarr(2*n_z+1) ;total density of non-thermal electrons (integrated over E and mu) at each z point, output
 S1=transpose(dblarr(n_E+1, n_mu+1, 2*n_z+1)) ;injection function (input), is assumed to be zero 
 
 t=0d0 
 
 t_arr=[t]
 f_arr=list(f)
 
 while t lt tmax do begin
  res=call_external(libname, 'SOLVE', t, S1, J, f, n_fast) ;perform one integration step
  print, t
  
  if (t-max(t_arr)) gt dt_record then begin
   t_arr=[t_arr, t]
   f_arr.add, f
  endif
 endwhile
 
 res=call_external(libname, 'DELETE_SOLVER', /unload) ;unload DLL
 
 ;compressing the data and converting them to GX Simulator-compatible format:
 N_E1=32 ;number of output energy nodes
 E=exp(alog(Emin)+(alog(Emax)-alog(Emin))*dindgen(N_E1)/(N_E1-1)) ;uniform log-spaced grid
 E_grid1=E_grid[0 : N_E-1] ;cut off the infinite energy point
 
 mu=mu_grid[0 : -1 : 2] ;select every 2nd pitch-angle node
 
 t=t_arr
 
 distfunc=dblarr(N_E1, n_elements(mu), n_elements(s), n_elements(t))
 
 dzx=(max(z_loop)-min(z_loop))/(n_elements(z_loop)-1) ;mean voxel size along the loop
 for ll=0, n_elements(t)-1 do begin
  f=transpose(f_arr[ll])
  f=f[0 : n_E-1, *, *] ;cut off the infinite energy point
  f=f[*, 0 : -1 : 2, *] ;select every 2nd pitch-angle node
  
  fx=dblarr(N_E, n_elements(mu), n_elements(s))
  for kk=0, n_elements(s)-1 do begin
   u=where(abs(z_grid-z_loop[kk]) lt (dzx/2))
   for ii=0, n_E-1 do for jj=0, n_elements(mu)-1 do fx[ii, jj, kk]=mean(f[ii, jj, u]) ;averaging the distribution function over z
  endfor
  
  for kk=0, n_elements(s)-1 do begin
   fxx=reform(fx[*, *, kk])
   u=where(fxx gt 0)
   fmin=min(fxx[u])
   u=where(fxx le 0, k)
   if k gt 0 then fxx[u]=fmin
   fx[*, *, kk]=fxx
  endfor
  
  for jj=0, n_elements(mu)-1 do for kk=0, n_elements(s)-1 do $
   distfunc[*, jj, kk, ll]=exp(interpol(alog(fx[*, jj, kk]), alog(E_grid1), alog(E)))
 endfor
 
 distfunc/=1d6 ;convert to CGS units
 s/=L ;convert the distance back to relative units
   
 save, distfunc, t, E, mu, s, L, filename='Evolution.sav', /compress
end