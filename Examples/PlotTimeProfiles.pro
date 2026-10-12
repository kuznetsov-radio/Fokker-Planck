pro PlotTimeProfiles
 restore, 'E:\Sci\Kinetic\Rozhkova\ExampleB.sav'
 
 N_E=n_elements(E)
 N_mu=n_elements(mu)
 N_t=n_elements(t)
 N_z=n_elements(z)
 
 s=z/L
 
 forward_function expnumber
 
 Elist=[0.01, 0.1, 1.0]
 
 slist=0.05+0.075*dindgen(13)
 
 nb=dblarr(n_elements(Elist), N_t, n_elements(slist))
 
 for k=0, n_elements(slist)-1 do for i=0, n_elements(Elist)-1 do begin
  a=min(abs(s-slist[k]), l)
  f1=reform(f[*, *, l, *])
  
  a=min(abs(E-Elist[i]), l)
  f1=reform(f1[l, *, *])
   
  for j=0, N_t-1 do nb[i, j, k]=2d0*!dpi*int_tabulated(mu, f1[*, j], /double)
 endfor
 
 set_plot, 'ps'
 
 for i=0, n_elements(Elist)-1 do begin
  device, file='TimeProfilesB'+string(i, format='(I02)')+'.eps', xsize=10.0, ysize=6.0, /encapsulated, /color, $
          bits_per_pixel=8, font_size=8
  !X.MARGIN=[7, 10]
  !Y.MARGIN=[3.2, 0.5]
  
  for k=0, n_elements(slist)-1 do begin
   nb_local=reform(nb[i, *, k])
   
   if k eq 0 then begin
    loadct, 0, /silent
    plot, t, nb_local, /ylog, yrange=max(nb[i, *, *], /nan)*[1.5d-4, 1.2], $
          ystyle=1, ytickformat='expnumber', $
          xtitle='!8t!3, s', ytitle='!8f!3', /nodata
    a=min(abs(E-Elist[i]), l)
    xyouts, 6400, 5500, '!8E=!3'+string(E[l], format='(F7.3)')+' MeV', /device
   endif else loadct, 34, /silent
   oplot, t, nb_local, color=255.0*k/(n_elements(slist)-1), thick=2
   xyouts, 8600, 5400-k*350, '!8z/L!3='+string(slist[k], format='(F5.3)'), /device, color=255.0*k/(n_elements(slist)-1)
  endfor   
  
  device, /close
 endfor 
end

function expnumber, axis, index, value
 u=floor(alog10(value))
 return, '10!U'+string(u, format='(2I0)')+'!N'
end