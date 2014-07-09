pro map_hps_ppt

  ; Program to take level p1-images (median subtracted), and produce p2-images, which are cleaned of HPs and CPs
  ;
  ; Notes: Only use images with have accurate JDates - work over the course of 1 orbit - 15 mins - gaps are ~ 85 mins.
  ;
  ; Produce mirror of _p1.sav in _p2.sav - with HPs and CPs taken care of.
  ; 
  ; Use maps
  ;
  Compile_opt idl2
  ;; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  
  ; input directory
  indir='~/BRITE/data/UB/testing/nostack/p1/'
  
  outdir='~/BRITE/data/UB/testing/nostack/p2/'
  
  filesin=file_search(indir+'*.sav', count=nsav)
  
  for kk=7, nsav-1 do begin
  
    ; restore observation file
    restore, filesin[kk]  ; roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp
    
    data2=data1[*,*,16:18]
    
    hdname=file_basename(filesin[kk], '_p1_test.sav')
    
    print, hdname
    
    ;restore hp map
    hpmap='~/BRITE/data/UB/reduction/hot_pixels/'+hdname+'_hpmap.sav'
    restore, hpmap  ; hp_xy
  
    npts=n_elements(jd) ; number of frames
    
    outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'
    fileout=outdir+hdname+'_hp_clean.ps'
    
    !p.multi=[0,2,1,0,0]
    ps_on, fileout, xsize=27, ysize=13, /landscape
      
      ; clean the image of HPs
      for j=0, 2, 2 do begin
        
        plotsym, 0, /fill, 0.4
       ; plot_image, bytscl(data2[*,*,j], 50, 500), color=cgcolor('black'), title=hdname, charsize=0.8
       ; oplot, hp_xy[0,*], hp_xy[1,*], psym=8, color=cgcolor('purple')
        
        frame0 = reform(data2[*,*,j])
        frame1 = B2_im2b(frame0,hp_xy)           ; image cleaned
        data2[*,*,j]=frame1
        
        ;window, 1, xsize=600, ysize=500, xpos=1100, ypos=300
        plot_image, bytscl(data2[*,*,j], 50, 500), color=cgcolor('black'), title=hdname+' - Time 1', charsize=0.8
        ;oplot, hp_xy[0,*], hp_xy[1,*], psym=2, color=cgcolor('orange')
        
;wait, 3
        
      endfor
      
      ps_off
      stop
      
     
      
    endfor
  
  
  stop
end
