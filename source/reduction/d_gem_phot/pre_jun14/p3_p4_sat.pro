pro p3_p4_sat

  ; modified from p3_p4_gem2 - to investigate saturated pixels
  
  ; 10 feb 2014
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, 0.8, /fill
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  indir='~/BRITE/data/UB/testing/p3_testsets/'
  ;window, 0, xsize=600, ysize=600
  
  hdname='37468'
  vmag='3.8'
    
  filesin=file_search(indir+'HD'+hdname+'*.sav', count=nfiles)
    
  fileout='~/BRITE/reports/myreports/deltadn_exptime_'+hdname+'.ps'
  ps_on, fileout, xsize=15, ysize=15
  plotsym, 0, /fill, 0.8
  cols=cgcolor(['forest green', 'purple', 'blue', 'orange red'])
  plot, [0,30], [0,2600.], color=cgcolor('black'), psym=8, xrange=[-1,31], title=hdname+', Vmag='+vmag, charsize=0.7, $
    xtitle='Consecutive pixels - from brightest', ytitle='Average DN in 1 orbit with error', /ynozero, /nodata
    
  for bb=0, 1 do begin  ; only do 0.1s and 1s ;nfiles-1 do begin
      
    print, filesin[bb]
      
    restore, filesin[bb] ;jd, data1
      
    jd1=jd-jd[0]
    if bb eq 0 then data1=data1*10.
    if bb eq 1 then data1=data1
    if bb eq 2 then data1=data1/3.
    if bb eq 3 then data1=data1/6.
      
    jd2=jd1[1:n_elements(jd1)-1]
    jdiff=jd2-jd1
      
    gap=where(jdiff gt 0.015, ngap)
    gap=[0,gap,n_elements(jd1)-1]
     
    ;    stop
      cadence=robust_mean(jdiff,2)  ; in days
      cadence=cadence*24.*60.       ; in minutes
      n_img_per_orbit=fix(15./cadence)
      
      nfrm=(size(data1, /dim))[2]
      
      med_pix=fltarr(30,ngap+1)
      scat_pix=fltarr(30,ngap+1)
      
      for gp=0, ngap do begin
        
        if gp eq 0 then iloc=indgen(gap[gp+1]+1) else iloc=indgen(gap[gp+1]-gap[gp])+gap[gp]+1  ; iloc is group of images in orbit
        
        nfrm=n_elements(iloc)
        
        data2=data1[*,*,iloc]
        
        dn1=fltarr(30,nfrm)
        
      for img=0, nfrm-1 do begin 
                
        dat1=reform(data2[*,*,img], 32*32)
         
        order1=reverse(sort(dat1))
        
        dn1[*,img]=(dat1[order1])[0:29]

      endfor   
; stop     
      for pix=0, 29 do med_pix[pix,gp]=median(dn1[pix,*])
      for pix=0, 29 do scat_pix[pix,gp]=robust_sigma(dn1[pix,*])
      
 ;    stop
     endfor
     
     med_pix1=fltarr(30)
     epix=fltarr(30)
     
     for i=0, 29 do med_pix1[i]=median(med_pix[i,*])
     for i=0, 29 do epix[i]=robust_mean(scat_pix[i,*],2)
    
       
     oplot, med_pix1, psym=8, color=cols[bb]
     oploterror, med_pix1, epix, psym=3, errcolor=cols[bb], hatlength=5, errthick=3
      
       
   endfor
   ;al_legend, ['0.1s', '1.0s', '3.0s', '6.0s'], textcolor=cgcolor('black'), psym=8, color=cols, /right, charsize=0.7
   al_legend, ['0.1s x 10', '1.0s'], textcolor=cgcolor('black'), psym=8, color=cols, /right, charsize=0.7
   
   ps_off
   
   spawn, 'open '+fileout+' &'
   stop
  

  stop
  print, 'end of program'
  stop
end


