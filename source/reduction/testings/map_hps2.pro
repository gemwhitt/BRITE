pro map_hps2, data1,exp_ttl, nhp,avhp

  ; Modified from map_hps
  ; 
  ; called by temp_trends
  ; 
  ; Purpose: find HPs in an individual frame - record number and average value 
  ; ; clean HPs then record average of top 5 target pixels
  ;
  ;
  ; CALLS: find_hps.pro
  ;
  Compile_opt idl2
  ;; FOR PLOTTING
  !p.background=cgcolor('white')
  
    ; number of frames in file
    nfrm=(size(data1, /dim))[2]
    
    ; loop over each image
    for ii=0, nfrm-1 do begin
    
      dat=data1[*,*,ii]/exp_ttl[ii]
      
      cr1=1.5 ;first criterion for rejection, e.g. 1.5*median to cope with increasing CCD temp...
      cr2=2   ;how many bad pixels
      
      ; check for HPs AND CPs
      ima=find_hps2(dat,cr1,cr2,ii, nhp,avhp)
      
    endfor
    
  print, 'End of map_hps2'
  
end