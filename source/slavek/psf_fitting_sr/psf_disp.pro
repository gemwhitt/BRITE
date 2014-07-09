pro psf_disp, im,tit,shx,shy,nbin
  ; SMR Oct.2013
  ;
  ; displays an 256x256 image
  ;
  ; uses the old version of TVimage.pro from Fanning
  ; requires : set_x
  ;
  ; input:
  ;   im = 256x256 pix image
  ;   tit = title, explanatory string
  ;   shx,shy = shifts in small pixels
  ;
  ; no output
  ; use:
  ;   B2_disp, psf0,'PSF-0',0,0
  ;
  ; may require erasing the display before use with: erase
  
 
  pos=[0.25,0.07,0.73,0.90]
  ; this requires adjustment for the position within the window
  
  plot,findgen(nbin*32.), $        ; needed to set hidden variables
    xsty=1,ysty=1, $
    pos=pos,/nodata
  erase                       ; stupid but must be
  
  TVimage2,bytscl(im),/keep,/over,pos=pos
  
  plot,findgen(nbin*32.), $
    xsty=9,ysty=9, $
    pos=pos,/nodata,/noerase
    
  oplot,[nbin*32./2.-10,nbin*32./2.+10]+shx,[nbin*32./2.,nbin*32./2.]+shy, color=cgcolor('orange')
  oplot,[nbin*32./2.,nbin*32./2.]+shx,[nbin*32./2.-10,nbin*32./2.+10]+shy, color=cgcolor('orange')
  
;  stop
  
  axis,xaxis=1,xran=!x.CRANGE/8,xsty=1,xtit=tit
  axis,yaxis=1,yran=!x.CRANGE/8,ysty=1
  
  return
end

