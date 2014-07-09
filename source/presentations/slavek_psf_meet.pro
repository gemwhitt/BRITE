pro slavek_psf_meet

  ; do psf-fitting for p3_test data - analyse the results
  ;
  Compile_opt idl2
  
  timefile='~/Desktop/timefile_2orbits.dat'
  ;
  ; FOR PLOTTING[[[[[[;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  ;
  todisplay='y'  ; see Step 5
  ;
  ;
  indir='~/BRITE/data/UB/testing/nostack/p3/'
  filesin=file_search(indir+'*.sav', count=nfiles)
  
    outdir='/Users/gemmawhittaker/BRITE/mini_meeting/plots/'

  norb=1  ; number of orbits (45 frames per orbit) to use to get average PSF
  
  for gg=0, nfiles-1 do begin
  
    print, file_basename(filesin[gg],'.sav')+' - file '+strtrim(gg,2)
;    stop
    
    restore, filesin[gg] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, sm_flx, avg_flx, expn
    
    fname=file_basename(filesin[gg],'.sav')
    dashpos=strpos(fname, '_')
    hdname=strmid(fname, 0, dashpos)
    
    dat=data1 ; no effect to data1
    
    ; get total number of frames in this file
    nfrm=(size(dat, /dim))[2]
    
    ; get number of orbits
    jd2=jd[1:nfrm-1]
    timediff=jd2-jd
    
    new_orbit=where(timediff gt 0.015, ngap)  ; there will be ngap+1 loops
    new_orbit=[new_orbit,nfrm-1]
    temp_orbit=new_orbit
    for rr=0+(norb-1), ngap, norb do temp_orbit[rr]=-999
    keep=where(temp_orbit eq -999)
    new_orbit=new_orbit[keep]
    ngap=n_elements(new_orbit)-1
    
    rbin=[8]
    
    for hh=0, 0 do begin
    
      ; skip step 1 (find HPs) and step 2 (clean HPs) - already done
      
      for orbit=-1, 0 do begin  ;ngap-1 do begin
      
        starttime=systime()
        
        if orbit eq -1 then dat1=dat[*,*,0:new_orbit[0]] else dat1=dat[*,*,new_orbit[orbit]+1:new_orbit[orbit+1]]
        
        ; rebin the images - 4x4, 8x8, 16x16...
        dat2=rebin_data(dat1,rbin[hh])
        
        ndat2=n_elements(dat2[0,0,*])
        
        ; find approx center of PSF
        mthr=1000
        spr=2*rbin[hh]
        x = fltarr(32*rbin[hh],32*rbin[hh])                   ; average image
        m = n_elements(dat2[0,0,*])          ; how many images
        for i=0,m-1 do x = x + dat2[*,*,i]   ; adding images
        x = x/m                             ; average image
        pks=brite_findpeaks(x, minthresh=mthr, spread=spr)
        
;        stop
        
        npks=n_elements(pks)/3.
        
        if npks gt 1 or pks[0] eq -999.9 then begin
          stop
          ; hot pixels remain - clean image again and add hp locations to map
          cr1=1000 ;first criterion for rejection, this is changed to 1.5*median to cope with increasing CCD temp...
          cr2=100 ;above local median
          cr3=2   ;how many bad pixels
          ima = B2_im1b(dat1,cr1,cr2,cr3, imc,ww)
        endif
        
        ; Step 3: Use a Gaussian as the 0-th approximation for the PSF...
        cen=1000.
        sig=2.
        psf0 = first_approx(cen,sig,rbin[hh],pks)
        
;        stop
;plotsym, 0, /fill, 0.8
;        fileout0=outdir+'slav_bin_pk.ps'
;        ps_on, fileout0, xsize=15, ysize=11
;        plot_image, bytscl(x, 50, 500), charsize=0.9, title=hdname
;        oplot, [pks[0]], [pks[1]], psym=8, color=cgcolor('orchid'), symsize=0.8
;        ps_off
        
;        fileout1=outdir+'slav_psf0.ps'
;        ps_on, fileout1, xsize=15, ysize=11
;        plot_image, bytscl(psf0, 50, 500), charsize=0.9, title=hdname
;        oplot, [pks[0]], [pks[1]], psym=8, color=cgcolor('orchid'), symsize=0.8
;        ps_off
        
 ;       stop
        
        
        
        ; Step 4: Determine the shifts in small pixels from the image centre for individual images assuming a given PSF.
        ; note that these shifts are given in "small pixels", so shx/8. is the shift in big pixels
        get_shifts, dat2,psf0,rbin[hh], shx1,shy1
        
        ; Step 5: Optional step - useful to see what is happening, e.g. one can check how well the centres are found.
        if todisplay eq 'y' then $
          for i=0, 0 do begin
          shx=shx1[i]
          shy=shy1[i]
          window, 0;, xpos=1500, ypos=200
          plot_image, dat2[*,*,i]
          
          window, 1;, xpos=2300, ypos=200
          plot_image, psf0
          stop
        endfor
        
        count=-1
        repeat begin
        
          count=count+1
          shx=shx1
          shy=shy1
          
          ; Step 6: Using the shx, shy, adjust individual images to the centre. The determine the improved 2-D PSF.
          ; Note the background is subtracted (for use as an PSF), but is given by "bck".
          psf1 = get_psf(dat2,shx,shy,rbin[hh], dud,bck,w)  ; w is 2d-indices of background pixels
          
          xw=(array_indices(findgen(256,256), w))[0,*]
          yw=(array_indices(findgen(256,256), w))[1,*]
          
          ;window, 2, xpos=2000, ypos=-300
          ;for i=0, n_elements(shx)-1 do begin
          ;plot_image, dat2[*,*,i]
          ;oplot, [shx[i]]+pks[0], [shy[i]]+pks[1], color=cgcolor('purple'), psym=2
          ;oplot, xw, yw, color=cgcolor('orange'), psym=1
          ;wait, 1
          ;endfor
          
          get_shifts, dat2,psf1,rbin[hh], shx1,shy1
          
          ; compute the residuals between shifts of old fit and shifts of new fit
          shx_res=abs(shx1-shx)
          shy_res=abs(shy1-shy)
          
          ; determine if there are any shifts and if so, repeat
          rs1=where(shx_res gt 0, nrs1)
          rs2=where(shy_res gt 0, nrs2)
          
          if nrs1 gt 0 OR nrs2 gt 0 then cycle=1 else cycle=0
          
        endrep until cycle eq 0  ; step 7 is to iterate until convergence
        
        ; count = number of iterations
        print, 'Number of iterations was '+strtrim(count,2)
        
        ; optional display....to see PSF
        if todisplay eq 'y' then begin
          psf_disp,psf1,'PSF-1',0.,0.,rbin[hh]
          stop
        endif
        
        ; Step 8: The flux & background determination by linear fits of the experimental PSF.
        psf3=psf1 ;[pks[0]-(5*rbin[hh]):pks[0]+(5*rbin[hh]),pks[1]-(5*rbin[hh]):pks[1]+(5*rbin[hh])]
        shx2=shx
        shy2=shy
        dud=dud
        yfit=1
        
        fileout2=outdir+'slav_psf1.ps'
                ps_on, fileout2, xsize=15, ysize=11
                plot_image, bytscl(psf1, 50, 500), charsize=0.9, title=hdname
                ps_off
               
        
        stop
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; this will be end of loop for finding the psf
        ;  stop
        fit = B2_im7b(psf3,dat2,shx2,shy2,dud, fite,yfit)
        ;stop
        for i=0, ndat2-1 do begin
        
          noise_map=(reform(dat2[*,*,0])*0.)+1.
          xloc=pks[0]+shx2[i]
          yloc=pks[1]+shy2[i]
          map_clean=1
          eflux=1
          mu=1
          
          window, 2, xpos=1500, ypos=-300
          plot_image, dat2[*,*,i]
          oplot, [xloc], [yloc], psym=2, color=cgcolor('orange')
          
          window, 3, xpos=2300, ypos=-300
          plot_image, psf3
          oplot, [xloc], [yloc], psym=2, color=cgcolor('orange')
          
          
          fastphot, reform(dat2[*,*,i]), psf3, noise_map, xloc, yloc, flux, map_clean=map_clean, eflux=eflux, mu=mu
          stop
          ; save fit and fite
          ;if orbit eq -1 then begin
          ;  flux=reform(fit[0,*])
          ;  ferr=reform(fite[0,*])
          ;  bkgd=reform(fit[1,*])
          ;  berr=reform(fite[1,*])
          ;endif else begin
          ;  flux=[flux,reform(fit[0,*])]
          ;  ferr=[ferr,reform(fite[0,*])]
          ;  bkgd=[bkgd,reform(fit[1,*])]
          ;  berr=[berr,reform(fite[1,*])]
          ;endelse
          ;stop
          ;imgsig=strmid(strtrim(robust_sigma(fit[0,*]),2), 0, 6)
          
          if orbit eq -1 AND i eq 0 then begin
            flux1=flux
            map_clean1=map_clean
          endif else begin
            flux1=[flux1,flux]
            map_clean1=[[map_clean1],[map_clean]]
          endelse
          stop
          endtime=systime()
        endfor
        
        ;openw, lun, timefile, /append, /get_lun
        ;printf, lun, hdname, rbin[hh], starttime, endtime, format='(a,2x,i,2x,a,2x,a)'
        ;free_lun, lun
        temp:
      endfor ; END LOOP OVER THIS SUN-SECTION OF ORBITS - IS USING
      
      ; save output
      ; name output file...
      ;outfile=outdir+hdname+'_rb'+strtrim(rbin[hh],2)+'_psf'+strtrim(norb,2)+'o.sav'
      ;save, filename=outfile, flux, ferr, bkgd, berr, data1, jd, roi_dim, xc, yc, ccd_temp, sm_flx, avg_flx, expn
      
    endfor  ; end loop over different bin sizes
    
    ;   stop
    
  endfor  ; end loop over file
  stop
end