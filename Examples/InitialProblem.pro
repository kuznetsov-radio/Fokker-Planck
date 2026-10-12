pro InitialProblem
 L=50000d3 ;loop length, m
 mr=5d0 ;mirror ratio
 
 n0=1d16 ;plasma density, m^{-3}
 T0=2.5d7 ;plasma temperature, K
 
 z0=25000d3/2 ;injection point, m
 w=2500d3 ;beam width, m

 nb0=1d8*1d6 ;initial electron density
 
 Emin=0.01d0 ;minimum energy, MeV
 Emax=2d0 ;maximum energy, MeV
 delta=2.5 ;power-law index
 
 ;model dimensions:
 n_E=60L
 n_mu=40L  ;must be even, the actual number of nodes will be n_mu+1
 n_z=1000L ;note that the actual number of nodes will be 2*n_z+1
 
 ;model switches:
 switch_turb=0L ;turbulence key
 switch_reverse_current=0L ;reverse current key
 
 ;model parameters:
 tmax=5d0        ;end time of the calculations, s
 dt0=0.01d0      ;initial integration time step (at t=0), s
 dt_max=0.01d0   ;maximum allowed integration time step, s
 dt_record=0.1d0 ;time interval between consequent saved states, s
 gamma=0.1d0     ;parameter that controls the integration accuracy and integration step (fractional accuracy)
 eps_z=0.25d0    ;z grid scaling factor
 
 libname='C:\pCloud\Codes\FokkerPlankSolver\x64\Release\FokkerPlankSolver.dll' ;Stupin's library name

;---------------------------------------------------
 
 ;fundamental constants (SI units):
 qe=1.602176634d-19
 me=9.1093837015d-31
 c=299792458d0
 eps0=8.8541878128d-12
 kB=1.380649d-23
 
 Nz0=1000
 z=L*dindgen(Nz0)/(Nz0-1) ;initial z array
 B=100d0*(1d0+(mr-1d0)*(2d0*z/L-1d0)^2) ;magnetic field profile
 
 ;plasma density and temperature profiles
 n_th=dblarr(Nz0)
 n_th[*]=n0
 T_th=dblarr(Nz0)
 T_th[*]=T0
 
 general_param=[n_E, n_mu, n_z, n_elements(z), switch_reverse_current, switch_turb] ;integer parameters of the model
 
 zmax=L ;maximum z coordinate
 r_gridz=eps_z*zmax ;parameter of the non-uniform z grid
 
 k_minE=1d0-exp(alog(1d0*n_E)/(1d0-Emax/Emin))
 r_gridE=-Emin/alog(1d0-k_minE) ;parameters of the non-uniform energy grid
 
 grid_param=[Emin, zmax, tmax, dt0, dt_max, gamma, r_gridE, r_gridz] ;floating-point parameters of the model
 
 coulomb_log=dblarr(n_elements(z)) ;Coulomb logarithm
 u=where(T_th lt 1.4d5, k)
 if k gt 0 then coulomb_log[u]=9.1-0.5*alog(n_th[u]/1d6)+1.5*alog(T_th[u])
 u=where(T_th ge 1.4d5, k)
 if k gt 0 then coulomb_log[u]=15.1-0.5*alog(n_th[u]/1d6)+alog(T_th[u])
 conductivity=3d0/4/sqrt(2d0*!dpi)*(4d0*!dpi*eps0)^2*(kB*T_th)^1.5/qe^2/sqrt(me)/coulomb_log ;conductivity [SI units]
 ;loop parameters:
 loop_param=transpose([transpose(z), transpose(B), transpose(n_th), transpose(coulomb_log), transpose(conductivity)]) 
 
 n_inf=1d20 ;chromospheric (outside the loop) plasma density, m^{-3}
 T_inf=10d3 ;chromospheric plasma temperature, K
 log_inf=(T_inf lt 1.4e5) ? 9.1-0.5*alog(n_inf/1d6)+1.5*alog(T_inf) : $
                            15.1-0.5*alog(n_inf/1d6)+alog(T_inf) ;Coulomb logarithm
 conductivity_inf=3d0/4/sqrt(2d0*!dpi)*(4d0*!dpi*eps0)^2*(kB*T_inf)^1.5/qe^2/sqrt(me)/log_inf ;conductivity [SI units]
 inf_param=[n_inf, log_inf, conductivity_inf] ;plasma parameters in the chromosphere
 
 turb_param=dblarr(n_E+1) ;mean free path - currently undefined
 turb_param[*]=1
  
 ;initialize the DLL using the above-mentioned parameters: 
 res=call_external(libname, 'CREATE_SOLVER', general_param, grid_param, inf_param, loop_param, turb_param)
 
 ;reserve space for the E, mu, z, etc. grids, create the irregular grids, and extrapolate the input data to these grids:
 E_grid=dblarr(n_E+1)
 mu_grid=dblarr(n_mu+1)
 z_grid=dblarr(2*n_z+1)
 dlnB_out=dblarr(2*n_z+1)
 n_loop_out=dblarr(2*n_z+1)
 res=call_external(libname, 'GET_GRIDS', E_grid, mu_grid, z_grid, dlnB_out, n_loop_out)  
 
 ;create the initial electron distribution function:
 f=dblarr(n_E+1, n_mu+1, 2*n_z+1) ;reserve space
 f_E=nb0*(delta-1)/Emin*(E_grid/Emin)^(-delta)
 f_mu=dblarr(n_mu+1)
 f_mu[*]=0.5
 f_z=exp(-(z_grid-z0)^2/w^2)
 for i=0, n_E do for j=0, n_mu do for k=0, 2*n_z do f[i, j, k]=f_E[i]*f_mu[j]*f_z[k]
 f=transpose(f) ;transpose the array for consistency with the DLL definitions          

 ;other input/output parameters: 
 Jr=dblarr(2*n_z+1) ;return current (input/output), is assumed to be zero
 n_fast=dblarr(2*n_z+1) ;total density of non-thermal electrons (integrated over E and mu) at each z point, output
 S=transpose(dblarr(n_E+1, n_mu+1, 2*n_z+1)) ;injection function (input), is assumed to be zero 
 
 t=0d0 
 
 t_arr=[t]
 f_arr=list(f)
 
 while t lt tmax do begin
  res=call_external(libname, 'SOLVE', t, S, Jr, f, n_fast) ;perform one integration step
  print, t
  
  if (t-max(t_arr)) gt dt_record then begin
   t_arr=[t_arr, t]
   f_arr.add, f
  endif
 endwhile
 
 res=call_external(libname, 'DELETE_SOLVER', /unload) ;unload DLL
  
 E=E_grid[0 : n_E-1] ;remove the highest (infinite) energy node
 mu=mu_grid
 
 uz=where((z_grid ge 0) and (z_grid le L))
 z=z_grid(uz) ;remove the spatial nodes outside the loop
 
 t=t_arr
 
 f=dblarr(n_E, n_mu+1, n_elements(z), n_elements(t))
 
 for i=0, n_elements(t)-1 do begin
  f1=transpose(f_arr[i])
  f1=f1[0 : n_E-1, *, *] ;remove the highest (infinite) energy node
  f1=f1[*, *, uz]        ;remove the spatial nodes outside the loop
  f[*, *, *, i]=f1
 endfor
 
 z=z[0 : -1 : 10]           ;save each 10th spatial node
 f=f[*, *, 0 : -1 : 10, *]  ;save each 10th spatial node
 
 save, f, E, mu, z, t, L, filename='E:\Sci\Kinetic\Rozhkova\ExampleA.sav', /compress
end