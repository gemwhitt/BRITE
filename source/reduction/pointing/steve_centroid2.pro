;+
;
; Name: steve_centroid2
;
; Purpose: Calculates the centroid of a star given an initial best
;          guess position
;
; Syntax: imstars2 = steve_centroid(imstars)
;
; Inputs: imstars - array of (x,y) positions of stars/features
;
; Opt Inputs: width - width of box to search for star. If width not
;                     set then default of 5 pixels is used.
;
; Outputs: imstars2 - array of better (x,y) positions of stars/features
;
; Calls: gcntrd
;
; Common: None
;
; Prev. Hist.: Algorithm originally developed by Steve Spreckley for
;              use with SMEI data
;
; History: 10-Oct-2007, Danielle Bewsher (RAL)
;          24-Oct-2007, Documented and made more elegant, DB
;
; Contact: Danielle Bewsher (d.bewsher@rl.ac.uk) and
;          Steve Spreckley (sas@star.sr.bham.ac.uk)
;
;-
FUNCTION steve_centroid2,data,imstars,width=width,fwhm=fwhm

  ;set box width to search for star peak
  IF (n_elements(width) eq 0) THEN width = 5. ; give a smaller width??
  
  dsz = size(data)
  
  x = round(imstars[0,*])
  y = round(imstars[1,*])
 
  isz = size(imstars,/dimensions)
  
  ;set up new array to collect results and output
  imstars2 = imstars*0.
  
  
  IF n_elements(isz) gt 1 THEN temp=isz[1] ELSE temp=1 ;so that even a single star can be tracked.
  
  FOR i=0,temp-1 DO BEGIN
    ;print,x,y
    
    IF (x[i]-width lt 0 or x[i]+width gt dsz[1]-1 or y[i]-width lt 0 or y[i]+width gt dsz[2]-1) THEN BEGIN
      imstars2[*,i] = imstars[*,i]
      
      stop
    ENDIF ELSE BEGIN
      
      sub = data[x[i]-width:x[i]+width,y[i]-width:y[i]+width]
      subx_guess = 16.
      suby_guess = 16.
      ;
      ;    plot_image,sub,origin=[x(i)-width,y(i)-width]
      ;    plots,x(i),y(i),psym=2,color=0
      
      bigsub = congrid(sub,33,33,/interp)
      
      bigsub2 = filter_image(bigsub,fwhm=fwhm,/all_pixels)  ; was 3.0,8 - Vino
      
     ; window, 1, xsize=500, ysize=500, xpos=1500, ypos=-500
     ; plot_image, bigsub
      
      ;wait, 3
      
      window, 2, xsize=500, ysize=500, xpos=2500, ypos=-500
      plot_image,bigsub2
      
      ;gcntrd,bigsub,subx_guess,suby_guess,xgcntrd3c,ygcntrd3c,fwhm,/silent   ; was 5.,8 -Vino
      centroid,bigsub,xgcntrd3c,ygcntrd3c,xy_peak=[subx_guess,suby_guess],fwhm=1.5
      
      
      ;wset, 1
      ;oplot, [xgcntrd3c], [ygcntrd3c],color=cgcolor('red'), psym=4, symsize=3
;      stop
      
      
      IF (xgcntrd3c ne -1 and ygcntrd3c ne -1) THEN BEGIN
        imstars2[0,i] = x[i]-width+xgcntrd3c/3.
        imstars2[1,i] = y[i]-width+ygcntrd3c/3.
        
        
        ;     plots,xgcntrd3c,ygcntrd3c,psym=1,color=1
        ;  wait,2.
      ENDIF ELSE BEGIN
        
        print, 'steve_centroid fail'
        
      ENDELSE
    ENDELSE
  ENDFOR
  
  return,imstars2
  
END
