pro plot_systematic

; investigate systematic trends - CCD temp etc

Compile_opt idl2
!p.background=cgcolor('white')
plotsym, 0, 0.4, /fill
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting

indir='~/BRITE/data/UB/p1/CENTAURUS/'  ; HP cleaned + p2_trend removed

outdir='~/BRITE/data/UB/plots/CENTAURUS/systematic_trends/'

filesin=file_search(indir+'*_p1b.sav', count=nsav)
fname=file_basename(filesin, '_p1b.sav')

for i=0, nsav-1 do begin
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                      ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                      
;                      stop
                      
  obj=obj_new('IDL_savefile', filesin[i])
  
  obj->restore, 'jd'
  obj->restore, 'data1'
  obj->restore, 'ccd_temp'
  obj->restore, 'medimg'
  obj->restore, 'roi_name'
  
  tname=roi_name[0]
  
jd=jd-jd[0]
nfrm=n_elements(jd)

; calculate median of each image after warm columns are removed
medimg2=fltarr(nfrm)
for j=0, nfrm-1 do medimg2[j]=median(data1[*,*,j])

jd2=jd[1:nfrm-1]
jdiff=jd2-jd
gap=where(jdiff gt 0.015, ngap)
gap=[-1,gap,nfrm-1]
ps_on, outdir+'medimg_raw/'+tname+'_medimgr.ps', xsize=16, ysize=15
for j=0, ngap do begin
  
  if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
  
  if j eq 0 then plot, medimg[iloc], color=cgcolor('black'), psym=8, xrange=[0,70], yrange=[50,350], ystyle=1, $
    xtitle='Exposure number in each orbit', ytitle='Image Median', title=tname+', raw', charsize=0.7 else $
    
  oplot, medimg[iloc], color=cgcolor('black'), psym=8

endfor
ps_off

ps_on, outdir+'medimg_p1b/'+tname+'_medimgp1.ps', xsize=16, ysize=15
for j=0, ngap do begin
  
  ;wset, 1

  if j eq 0 then iloc=indgen(gap[j+1]+1) else iloc=indgen(gap[j+1]-gap[j])+gap[j]+1
  
  if j eq 0 then plot, medimg2[iloc], color=cgcolor('black'), psym=8, xrange=[0.,70.], yrange=[-10.,90.], ystyle=1, $
    xtitle='Exposure number in each orbit', ytitle='Image Median', title=tname+', p1', charsize=0.7 else $
  
    oplot, medimg2[iloc], color=cgcolor('black'), psym=8
    
endfor
ps_off

endfor

print, 'end of program'
end


