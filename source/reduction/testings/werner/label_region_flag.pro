pro label_region_flag

; Program to determine frames with an adequate number of PSf pixels - use label_region.pro
; Use p2 images - which have been mostly - but not entirely cleaned of HPs
; 
; This modification of label_region_psf is for WERNER TESTSET files 
; 
; Output: Save a new variation of the input files, but in the date range of TOR and with flag added to denote
; good/bad/blank frames
;
Compile_opt idl2
  
indir='~/BRITE/TESTSETS/werner4lc/p2/'

filein=file_search(indir+'*CEN*129056*.sav', count=nf)

outdir='~/BRITE/TESTSETS/werner4lc/p2/subsets/'
  
for f=0, nf-1 do begin
  
  print, file_basename(filein[f], '.sav')
    
  restore, filein[f]  ;jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag
  
  ; shorten the testsets to the following date range (because this is the range for TOR)
  ; range is 27th June - 3rd July - 2014
  caldat, jd, mon, day, yr, hr, min, sec
  
  sday=27
  eday=3
  smon=6
  emon=7
  
  ; convert time range to JD
  jds=julday(smon,sday,2014)
  jde=julday(emon,eday,2014)
  
  iloc=where(jd ge jds AND jd le jde, ni)
     
  nimg=ni
  
  ; redefine arrays
  jd=jd[iloc]
  data1=data1[*,*,iloc]
  ccd_temp=ccd_temp[*,iloc]
  medcol1=medcol1[iloc]
  medcol2=medcol2[iloc]
  medimg0=medimg0[iloc]
  roi_loc=roi_loc[*,iloc]
  
  
  ; establish flag for ID
  flag=intarr(nimg)
  psf_loc=roi_loc*0
    
  ; set threshold for finding illuminated pixels
  thr=50
    
  nhp=intarr(nimg)
  npix_psf=intarr(nimg)
    
  for im=0, nimg-1 do begin
      
      ; find "bad" images and denote with flag=0
      if medimg0[im] le 20 OR medimg0[im] ge 5000 then begin
        print, 'Bad image '
       
        ;plot_image, bytscl(data1[*,*,im], 20, 200)
        
        goto, skipimage
      endif
 
      dat=data1[*,*,im]
      
      ; make all pixels in dat=1 equal to 0 instead
      xx=where(dat eq 1, nxx)
      if nxx gt 0 then begin
        twod=array_indices(dat, xx)
        dat[twod[0,*],twod[1,*]]=0
      endif
      
      s=size(dat, /dim)
      
      ; add a border to the image
      dat2=lonarr(s[0]+2,s[1]+2)
      dat2[1,1]=dat
      
      ; use label region to determine number of illuminated regions and number of PSF pixels
      r2=label_region(dat2 ge thr, /all_neighbors)
      r2=r2[1:s[0],1:s[1]]                         ; trim off border
      
      hist2=histogram(r2, binsize=1, locations=loc2, reverse_indices=ri2)
      
      if n_elements(hist2) eq 1 then begin
        flag[im]=1
        goto, skipimage
      endif
      
;xx=where(hist2 eq 1, nxx)
;plot_image, bytscl(dat, 20, 200)
;for ii=0, nxx-1 do begin
;  ind=ri2[ri2[xx[ii]]:ri2[xx[ii]+1]-1]
;  ind2=array_indices(dat, ind)
;  oplot, [ind2[0,*]], [ind2[1,*]], color=cgcolor('green'), psym=2
;endfor
      
      ; the number of results in hist2 are the number of regions (blobs) > thr, plus the background pixels
      ; so if n_elements(hist2) > 2 then - HP!
      ; record n_elements(hist2)-1 + number of pixels in second largest "blob" - i.e. PSF pixels
      nhp[im]=n_elements(hist2)-2
      

      ; now overplot 2nd biggest group of pixels
      ; sort hist1 in decreasing order
      sort2=reverse(sort(hist2))
      ind2=ri2[ri2[sort2[1]]:ri2[sort2[1]+1]-1]
      i2=array_indices(dat, ind2)
      
      xi=i2[0,*]
      yi=i2[1,*]
      
      ; get borders of PSF
      psf_loc[*,im]=[min(xi),max(xi),min(yi),max(yi)]
;wset, 0
;
;oplot, i2[0,*], i2[1,*], color=cgcolor('green'), psym=2
      
      npix_psf[im]=hist2[sort2[1]]
      ;stop
      skipimage:
    endfor  ; end loop over image
    
    ; remove zeros
    n2=npix_psf[where(npix_psf gt 5)]
    
    plot, n2, color=cgcolor('black'), psym=2
    res=n2-smooth(n2, 45, /edge_truncate)
    plot, res, color=cgcolor('purple'), psym=2
    
    xx=where(res ge median(res)-(robust_sigma(res)*3) AND $
      res le median(res)+(robust_sigma(res)*3), nxx)
    
    ; images at xx are "good PSFs" - flag=2  
    flag[xx]=2
    ; save output
    fileout=outdir+file_basename(filein[f])
    
    save, filename=fileout, jd, data1, ccd_temp, medcol1, medcol2, medimg0, roi_loc, vmag, bmag, flag, psf_loc
    
  endfor  ; end loop over file
  
  
  print, 'end of program'
end