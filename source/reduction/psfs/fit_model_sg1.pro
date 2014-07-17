pro fit_model_sg1

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
  
  for im=0, nfrm-1 do begin
    
    xdim=(size(model2, /dim))[0]
    ydim=(size(model2, /dim))[1]
    
    data=reform(data1[*,*,im])
    
    ; find where max pixel in image is
    imax=where(data eq max(data))
    xy_max=array_indices(data, imax)
    
    ; do data2-model - shift the model using xy_psf
    mod2=model2*0.0
    xymod=where(model2 gt 0)
    xmod=(array_indices(model2, xymod))[0,*]
    ymod=(array_indices(model2, xymod))[1,*]
    ;
    xc=xdim/2.
    yc=ydim/2.
    xdif=xc-xy_max[0]
    ydif=yc-xy_max[1]
    ;
    mod2[xmod-xdif,ymod-ydif]=model2[xmod,ymod]
    
    data=float(data)
    
    data2=data/max(data)
    
    res=data2-mod2
    
    wset, 0
    plot_image, bytscl(res, -1, 1)
    
    stop
    
    
    
    
  endfor
  
  
  
  
  
  stop
  
  
endfor



print, 'end of program'
end