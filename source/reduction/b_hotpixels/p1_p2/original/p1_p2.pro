pro p1_p2

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
toplot='n'

nstk=5
  
; input directory
indir='~/BRITE/data/UB/p1/CENTAURUS/'+strtrim(nstk,2)+'stk/'
  
outdir='~/BRITE/data/UB/p2/CENTAURUS/'+strtrim(nstk,2)+'stk/'

hpdir='~/BRITE/data/UB/reduction/hot_pixel_maps/CENTAURUS/'+strtrim(nstk,2)+'stk/'
  
filesin=file_search(indir+'*.sav', count=nsav)
fname=file_basename(filesin, '.sav')
  
for kk=0, nsav-1 do begin
  
  ; restore observation file
  restore, filesin[kk]  ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                        ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                        
                        ;jd, jd1, data1, ccd_temp, roi_dim, simbad_mags - testsets
                        
  print, file_basename(filesin[kk],'.sav')
  
  hdpos=strpos(fname[kk], 'HD')
  hdname=strmid(fname[kk], hdpos)

  ;restore hp map
  hpmap=hpdir+fname[kk]+'_hp_map.sav'
  
  restore, hpmap  ; hp_xy, roi_dim
  hp_xy=hp_xy[*,where(hp_xy[0,*] ne -1)]

  nfrm=n_elements(jd) ; number of frames
  
  jd1=jd-jd[0]
    
  ; clean the image of HPs
  for j=0, nfrm-1 do begin
      
    img0 = reform(data1[*,*,j])
    img1 = clean_img(img0,hp_xy)           ; image cleaned of HPs and CPs
    data1[*,*,j]=img1
    
    if toplot eq 'y' then begin
    if j eq 0 OR j eq nfrm-1 then begin
      window, 0, xsize=600, ysize=500, xpos=200, ypos=300
      plot_image, bytscl(img0, 50, 500)
      oplot, hp_xy[0,*], hp_xy[1,*], psym=2, color=cgcolor('orange')
      
      window, 1, xsize=600, ysize=500, xpos=800, ypos=300
      plot_image, bytscl(img1, 50, 500)
     
      stop
    endif
    endif  
  endfor  ; end loop over images
  
;  exp_time=exp_time[0]
;  exp_ttl=exp_ttl[0]
;  stop
  ; save new file
  fn=file_basename(fname[kk],'_p1_'+strtrim(nstk,2))
  fileout=outdir+fn+'_p2_'+strtrim(nstk,2)+'.sav'
  
  save, filename=fileout, roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                        medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
    
     
endfor  ; end loop over file
  
print, 'end of program'
end
