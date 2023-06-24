pro PlotLC
 restore, 'Evolution_noTurbulence.sav'
 
 N_E=n_elements(E)
 N_mu=n_elements(mu)
 N_t=n_elements(t)
 N_s=n_elements(s)
 
 s_list=[-0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
 N_loc=n_elements(s_list)
 
 E1list=[0.01,  0.01, 0.025, 0.05,  0.1, 0.25, 0.5, 1.0]
 E2list=[ 2.0, 0.025,  0.05,  0.1, 0.25,  0.5, 1.0, 2.0]
 c_list=[1.0, 1.0, 15.0, 100.0, 800.0, 10000.0, 80000.0, 6d5]
 N_int=n_elements(E1list)
 
 forward_function nb_total
 
 nb=dblarr(N_t, N_s, N_int)
 
 for i=0, N_t-1 do for j=0, N_s-1 do begin
  f=reform(distfunc[*, *, j, i])
  
  for k=0, N_int-1 do begin
   u=where((E ge E1list[k]) and (E le E2list[k]))
   f1=f[u, *]
   E1=E[u]
   nb[i, j, k]=nb_total(f1, E1, mu)
  endfor 
 endfor
 
 set_plot, 'ps'
 
 for k=0, N_int-1 do begin
  nbl=reform(nb[*, *, k])
  
  device, file='DS'+string(k, format='(I02)')+'.eps', xsize=8.9, ysize=6, font_size=8, bits_per_pixel=8, $
          /encapsulated, /color
  !X.MARGIN=[8, 1]
  !Y.MARGIN=[3.2, 0.5]
  loadct, 3, /silent
  contour, alog10(nbl), t, s, xstyle=1, ystyle=1, $
           xtitle='!17Time!3, s', ytitle='!17Distance, relative to the loop length!3', $
           levels=max(alog10(nbl))+[-1000, 5d0*(dindgen(100)/100-1)], /fill
  contour, alog10(nbl), t, s, levels=max(alog10(nbl))+[-5, -3, -2, -1], /overplot, c_color=255   
  xyouts, 100, 100, '!6'+string(E1list[k], format='(F5.3)')+' - '+string(E2list[k], format='(F5.3)')+' MeV!3', /device
  device, /close
 endfor     
 
 for l=0, N_loc-1 do begin
  a=min(abs(s-s_list[l]), j)
  nbl=reform(nb[*, j, *])
  
  ymax=-1d100
  ymin=1d100
  mf=0
  for k=0, N_int-1 do begin
   lc=reform(nbl[*, k])*c_list[k]
   
   ymax=ymax>max(lc)
   
   y=min(lc, xc)
   if xc eq 0 then begin
    y=lc[N_t-1]
    mf=1
   endif
   
   ymin=ymin<y
  endfor
  ymax*=2
  ymin/=(mf ? 10 : 2)
  
  device, file='LC'+string(l, format='(I02)')+'.eps', xsize=8.9, ysize=6, font_size=8, bits_per_pixel=8, $
          /encapsulated, /color
  !X.MARGIN=[8, 1]
  !Y.MARGIN=[3.2, 0.5]
  loadct, 39, /silent
  for k=0, N_int-1 do begin
   lc=reform(nbl[*, k])*c_list[k]
   if k eq 0 then plot, t, lc, /nodata, xstyle=1, /ylog, ystyle=1, yrange=[ymin, ymax], $
                        xtitle='!17Time!3, s', ytitle='!17Density!3, cm!U-3!N'
   oplot, t, lc, color=255.0*k/(N_int-0.9), thick=(k eq 0) ? 2 : 1    
   xyouts, 1800, 2300-k*200, string(E1list[k], format='(F5.3)')+'-'+string(E2list[k], format='(F5.3)')+$
           ' (!9X!3'+string(c_list[k], format='(I0)')+')', $
           /device, color=255.0*k/(N_int-0.9), charsize=6.0/8
  endfor
  xyouts, 100, 100, '!18s!6='+string(s_list[l], format='(F5.2)')+'!3', /device
  device, /close
 endfor    
end

function nb_total, f, E, mu
 N_E=n_elements(E)
 
 Q=dblarr(N_E)
 for i=0, N_E-1 do Q[i]=2d0*!dpi*int_tabulated(mu, f[i, *], /double)
 
 return, int_tabulated(alog(E), Q*E, /double)
end
