pro thr_aper, model, nbin, thr_aper1, flag

  plotsym, 0, /fill, 0.5
  
  ; program to take model PSF and use the label region with a threshols to determine best aperture
  ;
  ; Called by psf_aper
  ;
  ; Output: psf_aper[*,*,3] fo r3 apertures where *,* are dimensions of ccd raster
  ;
  ;
  ; Old method.....
        thr=[100,50,25,12]
        
        thr_aper1=[]
        
        for th=0, n_elements(thr)-1 do begin
          ;
          ;      ; use label region to determine number of illuminated regions and number of PSF pixels
          r1=label_region(model ge thr[th], /all_neighbors)
          
          hist1=histogram(r1, binsize=1, locations=loc1, reverse_indices=ri1)
          
          if n_elements(hist1) eq 1 then continue 
          
          sb=size(model, /dim)
          
          ; sort hist2 in decreasing order
          sort1=reverse(sort(hist1))
          
          ind=ri1[ri1[sort1[1]]:ri1[sort1[1]+1]-1]
          i1=array_indices(model, ind)
          
          xi=reform(i1[0,*])
          yi=reform(i1[1,*])
          
;            plot_image, bytscl(model, 0, 200)
;                  oplot, xi, yi, color=cgcolor('purple'), psym=2
;                  stop
;          
          thr_aper=fltarr(sb[0],sb[1])
          thr_aper[xi,yi]=model[xi,yi]
          
          thr_aper1=[[[thr_aper1]],[[thr_aper]]]

        endfor
 
  ;
  ;
end