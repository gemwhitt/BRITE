pro map_hps

; Modified from map_hps_sep
; Purpose: For each ROI save-file, make a map of the HP locations only....
; ... use the whole sequence to find HPs which are sometimes obscured by the target flux
; 
; 
; CALLS: find_hps.pro
;
Compile_opt idl2
;; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')
plotsym, 0, /fill, 1.

; input directory
indir='~/BRITE/data/UB/p1/ORION/'
  
outdir='~/BRITE/data/UB/reduction/hot_pixel_maps/ORION/'
  
filesin=file_search(indir+'*.sav', count=nsav)
fname=file_basename(filesin, '.sav')

if nsav eq 0 then stop
  
for bb=0, nsav-1 do begin
  
  hdname=file_basename(filesin[bb], '_p1.sav') 
  
  print, hdname 
  
  ; restore data
  obj=obj_new('IDL_Savefile', filesin[bb])
  obj->restore, 'data1'
  obj->restore, 'jd'
  obj->restore, 'medimg0'
  obj->restore, 'ccd_temp'
  
  nfrm=n_elements(jd)
  
  totdn=lonarr(nfrm)
  for img=0, nfrm-1 do totdn[img]=total(data1[*,*,img])
  
  rej1=where(medimg0 gt 5000 OR totdn lt 5000, nrej1, complement=keep1) ; bad images - discard
  
  nkeep=n_elements(keep1)
  
  jd=jd[keep1]
  data1=data1[*,*,keep1]
  medimg0=medimg0[keep1]
  totdn=totdn[keep1]
  ccd_temp=ccd_temp[*,keep1]
  
  nfrm=n_elements(jd)
  
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
  
  dat1=data1
   
  ; check to see if hpmap exists for this ROI
  outfile=outdir+fname[bb]+'_hpmap.sav'
  chk_file=file_search(outfile, count=chk)
  
  if chk eq 1 then restore, outfile  ; gives hp_xy
  
  ; determine number of orbits - 1 orbit=15 mins (0.01 days), 1 gap = 85 mins (0.059 days)
  jd1=jd-jd[0]
  jd2=jd1[1:nfrm-1]
  jdiff=jd2-jd1
  gap=where(jdiff gt 0.015, ngap)  ; there are ngap+1 observations
  gap=[-1,gap,nfrm-1]
  
  num_hp=lonarr(ngap+1)
  avg_temp=fltarr(ngap+1)
    
  ; loop over each observation window
  for ii=0, ngap-1 do begin
    
    iloc=indgen(gap[ii+1]-gap[ii])+gap[ii]+1
    
    dat=dat1[*,*,iloc]
    npts=n_elements(iloc)
    
    if npts lt 2 then continue

    cr1=50.+median(dat)   ;cr1 is proportional to the median of the observations
    cr2=2   ;how many bad pixels

    ; check for HPs AND CPs
    ima=find_hps(dat,cr1,cr2, ww,wx,wy)
    
    nhp=n_elements(ww)
    num_hp[ii]=nhp
    avg_temp[ii]=avg(ccd_temp[*,iloc])
    
    ;plot_image, bytscl(ima, 20, 100)
    ;oplot, [wx], [wy], psym=8, color=cgcolor('purple')
    ;wait, 0.4
    
    if nhp eq 0 then goto, no_hp
;    stop
    hp_loc=intarr(2,nhp)
    hp_loc[0,*]=wx
    hp_loc[1,*]=wy   
    
    ; check if there are new HPs and if so - save them to the map
    if n_elements(hp_xy) eq 0 then hp_xy=hp_loc else hp_xy=[[hp_xy],[hp_loc]]
    
    ;wset, 0
    ;plot_image, bytscl(ima, 20, 300)
    ;oplot, hp_xy[0,*], hp_xy[1,*], color=cgcolor('orange'), psym=2
    
;stop
    
    
    no_hp:
  endfor  ; end loop over each orbit
  
  
  plot, avg_temp, num_hp, color=cgcolor('purple'), psym=8
  stop
  ; remove duplicates
  if n_elements(hp_xy) gt 0 then begin
    
  all_xy=reform(strmid(strtrim(hp_xy[0,*],2),0)+'_'+strmid(strtrim(hp_xy[1,*],2),0))
            
  uloc=uniq(all_xy, sort(all_xy))
          
  hp_xy=hp_xy[*,uloc]
  
  save, filename=outfile, hp_xy, roi_dim
  
  undefine, hp_xy
  
  endif 
  
endfor  ; end loop over file (roi)

print, 'End of program'

end