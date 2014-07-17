pro map_hps_orion

  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  
  ; test method to correct HPs in individual images
  ;
  Compile_opt idl2
  !p.background=cgcolor('white')
  plotsym, 0, /fill, 0.9
  ;window, 0, xsize=600, ysize=550, xpos=1000, ypos=200
  
  ; use test data
  indir='~/BRITE/data/UB/p1/ORION/'
  
  outdir='~/BRITE/data/UB/reduction/hot_pixel_maps/ORION/'
  
  filesin=file_search(indir+'*.sav', count=nf)
  
  for f=0, nf-1 do begin
  
    restore, filesin[f] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, ccd_temp, exp_time, exp_ttl, $
    ;medcols, medimg0, medimg1, ndead, nsat, simbad_radec, vmag, bmag, parlax, otype, sptype
    
    nfrm=n_elements(jd)
    
    fname=file_basename(filesin[f],'_p1.sav')
    
    print, fname
    print, vmag
    
    totdn=lonarr(nfrm)
    for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
    
    rej1=where(medimg0 gt 5000 OR totdn lt 5000, nrej1, complement=keep1) ; bad images - discard
    
    nkeep=n_elements(keep1)
    
    jd=jd[keep1]
    data1=data1[*,*,keep1]
    ccd_temp=ccd_temp[*,keep1]
    medimg0=medimg0[keep1]
    totdn=totdn[keep1]
    
    xdim=(size(data1, /dim))[0]
    ydim=(size(data1, /dim))[1]
    
    hp1=lonarr(xdim, ydim, nkeep)
    
    jd1=jd-jd[0]
    jd2=jd1[1:nkeep-1]
    jdiff=jd2-jd1
    gap=where(jdiff gt 0.015, ngap) ; ngap+1 orbits
    gap=[-1,gap,nkeep-1]
    
    for orbit=5, ngap do begin
    
      if orbit mod 100 eq 0 then print, orbit
      
      
      iloc=indgen(gap[orbit+1]-gap[orbit])+gap[orbit]+1
      
      stop
      
      dat=data1[*,*,iloc]
      nimg=n_elements(iloc)
      
      thr=[]
      
      for img=0, nimg-1 do begin
      
        data2=dat[*,*,img]
        
        dead=where(data2 eq 0.0 OR data2 ge 8000, ndead, complement=live)
        
        sat=where(data2 ge 8000, nsatp)
        if nsatp gt 0 then begin
          xs=(array_indices(data2, sat))[0,*]
          ys=(array_indices(data2, sat))[1,*]
          data2[xs,ys]=7999
        endif
        
        nlive=n_elements(live)
        
        xl=(array_indices(data2, live))[0,*]
        yl=(array_indices(data2, live))[1,*]
        
        ; get median and scatter of live pixels
        meanl=mean(data2[xl,yl])
        scatl=robust_sigma(data2[xl,yl])
        
        cr1=(meanl+3*scatl < 100)
        
        cr2=4
        if cr1 le 0.0 then stop
        
        thr=[thr,cr1]
        
        ima=find_hps2(data2,cr1,cr2,iloc[img], ww,wx,wy)
        
        ; plot_image, bytscl(data2, 20, 200)
        ; oplot, [wx], [wy], color=cgcolor('orange'), psym=2
        
        nnew=n_elements(ww)
        
        if nnew eq 0 then goto, next_image
        
        data2=data2
        
        ; record HP locations and values in each image
        if n_elements(hps) eq 0 then begin
        
          nhp=nnew
          
          hps=fltarr(3,nhp,nimg)
          
          for hp=0, nhp-1 do begin
            hps[0, hp, *]=wx[hp]
            hps[1, hp, *]=wy[hp]
            hps[2, hp, img]=data2[wx[hp],wy[hp]]
          endfor
          
          
        endif else begin
        
          for hp=0, nnew-1 do begin
          
            ; check if HP found previously
            
            hp_match=where(hps[0,*,0] eq wx[hp] AND hps[1,*,0] eq wy[hp], nmatch)
            
            if nmatch gt 1 then stop
            
            if nmatch eq 1 then begin
            
              hps[2,hp_match, img]=data2[wx[hp],wy[hp]]
              
            endif
            
            if nmatch eq 0 then begin ; new entry - grow array
            
              nhp=(size(hps, /dim))[1]
              
              hp_temp=fltarr(3,nhp+1,nimg)
              
              for oldhp=0, nhp-1 do begin
                hp_temp[0, oldhp, *]=hps[0,oldhp,*]
                hp_temp[1, oldhp, *]=hps[1,oldhp,*]
                hp_temp[2, oldhp, *]=hps[2,oldhp,*]
              endfor
              
              hp_temp[0,nhp,*]=wx[hp]
              hp_temp[1,nhp,*]=wy[hp]
              hp_temp[2,nhp,img]=data2[wx[hp],wy[hp]]
              
              hps=hp_temp
              
            endif
            
          endfor
          
        endelse
        
        next_image:
      endfor  ; end loop over this image
      
      if n_elements(hps) eq 0 then goto, next_orbit
      
      ; check frequency of HP occurance
      nhp=(size(hps, /dim))[1]
      count=nhp
      
      if nimg lt 2 then begin
      
        hps=reform(hps, 3,nhp,1)
        hps2=reform(hps, 3,nhp,1)
      endif else hps2=hps
      
      
      for hp=0, nhp-1 do begin
      
        xx=where(hps[2,hp,*] gt thr, nxx)
        
        if hps[0,hp,0] eq 12 AND hps[1,hp,0] eq 0 then stop
        
        if nxx lt round(0.5*nimg) then begin  ; not hp
          count=count-1
          hps2[0:2,hp,*]=-999.9  ; get rid of hp
        endif
       
      endfor
      stop
      if count lt 1 then goto, next_orbit
      
      newhps=fltarr(3,count,nimg)
      
      newhps[0,*,*]=hps[0,where(hps2[0,*,0] ne -999.9),*]
      newhps[1,*,*]=hps[1,where(hps2[0,*,0] ne -999.9),*]
      newhps[2,*,*]=hps[2,where(hps2[0,*,0] ne -999.9),*]
      
      hps=newhps
      
      ;  for kk=0, nimg-1 do begin
      ;    plot_image, bytscl(dat[*,*,kk], 20, 200)
      ;    oplot, hps[0,*,kk], hps[1,*,kk], psym=2, color=cgcolor('blue')
      ;
      ;  wait, 0.5
      ;  endfor
      
      ; save HP locations for each image - to a .sav file
      ; update arrays
      wx1=fix(reform(hps[0,*,0]))
      wy1=fix(reform(hps[1,*,0]))
      
      for img=0, nimg-1 do begin
      
        temp1=lonarr(xdim,ydim)
        temp1[wx1,wy1]=1
        
        hp1[*,*,iloc[img]]=temp1
        
      endfor
      
      next_orbit:
      undefine, hps
      
    endfor ; end loop over this orbit
    
   ; plotsym, 0, /fill, 1.3
   ; for i=1000, 1100 do begin
   ;   plot_image, bytscl(data1[*,*,i], 20, 200)
   ;   hh=where(hp1[*,*,i] eq 1)
   ;  hx=(array_indices(hp1[*,*,i], hh))[0,*]
   ;  hy=(array_indices(hp1[*,*,i], hh))[1,*]
   ;   oplot, hx,hy, color=cgcolor('blue'), psym=8
   ;   wait, 0.2
   ; endfor
    
    fileout=outdir+fname+'_hpmap.sav'
    print, fileout
    
    save, filename=fileout, hp1
    
  endfor  ;end loop over this file
  
  
  
  print, 'end of program'
end