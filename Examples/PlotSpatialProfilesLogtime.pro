pro PlotSpatialProfilesLogtime
 restore, 'E:\Sci\Kinetic\Rozhkova\ExampleA.sav'
 
 N_E=n_elements(E)
 N_mu=n_elements(mu)
 N_t=n_elements(t)
 N_z=n_elements(z)
 
 s=z/L
 
 forward_function expnumber
 
 Elist=[0.01, 0.1, 1.0]
 
 tlist=[0.0, exp(alog(0.1)+(alog(5.0)-alog(0.1))*dindgen(12)/11)]
 
 nb=dblarr(n_elements(Elist), N_z, n_elements(tlist))

 for k=0, n_elements(tlist)-1 do for i=0, n_elements(Elist)-1 do begin
  a=min(abs(t-tlist[k]), l)
  f1=reform(f[*, *, *, l])

  a=min(abs(E-Elist[i]), l)
  f1=reform(f1[l, *, *])

  for j=0, N_z-1 do nb[i, j, k]=2d0*!dpi*int_tabulated(mu, f1[*, j], /double)
 endfor
 
 set_plot, 'ps'
 
 for i=0, n_elements(Elist)-1 do begin
  device, file='SpatialProfilesA'+string(i, format='(I02)')+'.eps', xsize=10.0, ysize=6.0, /encapsulated, /color, $
          bits_per_pixel=8, font_size=8
  !X.MARGIN=[7, 10]
  !Y.MARGIN=[3.2, 0.5]  
  
  for j=0, n_elements(tlist)-1 do begin
   nb_local=reform(nb[i, *, j])
    
   if j eq 0 then begin
    loadct, 0, /silent
    plot, s, nb_local, xstyle=1, /ylog, yrange=max(nb[i, *, *], /nan)*[1.5d-4, 1.2], $
          ystyle=1, ytickformat='expnumber', $
          xtitle='!8z!3/!8L!3', ytitle='!8f!3', /nodata
    a=min(abs(E-Elist[i]), l)
    xyouts, 6400, 5500, '!8E=!3'+string(E[l], format='(F7.3)')+' MeV', /device
   endif else loadct, 34, /silent
   
   oplot, s, nb_local, color=255.0*j/(n_elements(tlist)-1), thick=2
   xyouts, 8800, 5400-j*350, '!8t!3='+string(tlist[j], format='(F4.2)')+' s', /device, color=255.0*j/(n_elements(tlist)-1)
  endfor
  
  device, /close
 endfor
end

function expnumber, axis, index, value
 u=floor(alog10(value))
 return, '10!U'+string(u, format='(2I0)')+'!N'
end