pro estimate_flux

  ; fit models to sg1 data - get stats
  
  Compile_opt idl2
  !p.background=cgcolor('white')
  
  indir='~/BRITE/data/UB/p2/ORION/sg1/'
  
  fin=file_search(indir+'sav/*.sav', count=nf)
  
  if nf eq 0 then stop
  
  for ff=0, nf-1 do begin
  
    restore, fin[ff]  ;vmag, hdname, jd, data1, ccd_temp, xy_psf
    
    restore, indir+'models/'+hdname+'_mod1.sav' ; model, model2 (normalised)
    
    nfrm=n_elements(jd)
    
    value=fltarr(nfrm)
    res=fltarr(nfrm)
    
    for im=0, nfrm-1 do begin
    
     
      data=reform(data1[*,*,im])
      
      ; find where max pixel in image is
      imax=where(data eq max(data))
      xy_max=array_indices(data, imax)  ; 2-element array
      
      ; find where model2 gt 0
      mod2=model2*0.0
      xymod=where(model2 gt 0)
      xmod=(array_indices(model2, xymod))[0,*]
      ymod=(array_indices(model2, xymod))[1,*]
      
      ; find where max model is
      max_mod=array_indices(model, where(model eq max(model))) ; 2-element array
      
      ; calcualte difference between model and image peaks
      difx=max_mod[0]-xy_max[0]
      dify=max_mod[1]-xy_max[1]

      ; get total_flx in model
      flx_mod=total(model)
      
      ; get total flx in image based on location of model pixels
      flx_img=total(data[xmod-difx,ymod-dify])
      
      ;wset, 0
      ;plot_image, bytscl(data, 20, 500)
      ;oplot, [xmod-difx], [ymod-dify], color=cgcolor('orange'), psym=2
      
      flx_res=abs(flx_img-flx_mod)
      
      pc_of_im=flx_res/flx_mod*100.
      
      value[im]=pc_of_im
      
      res[im]=flx_res
     
    endfor
    
    wset, 0
    plot, jd, res, color=cgcolor('black'), psym=2
    
    wset, 1
    plot, jd, value, color=cgcolor('black'), psym=2

    
    
    stop
    
    
  endfor
  
  
  
  print, 'end of program'
end