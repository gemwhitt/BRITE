pro p1_info

; get data characteristics from p1 data for each month, e.g. number of GOOD observations versis BAD...
; CCD_TEMP changes, avg_median, HPs
; 
Compile_opt idl2
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')

indir='/Users/gemmawhittaker/BRITE/data/UB/p1/ORION/'

outdir='/Users/gemmawhittaker/BRITE/data/UB/reduction/p1_info/'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

fname=file_basename(filesin, '_p1.sav')

dashpos=strsplit(fname, '_')

hdname=strarr(nf)
for i=0, nf-1 do hdname[i]=strmid(fname[i], dashpos[i,1])

; get uniq hdnames
uname=hdname[uniq(hdname, sort(hdname))]
nuniq=n_elements(uname)

alltemps=[]
alljds=[]
allttls=[]
imgmeds=[]
imgmeds2=[]

st=dblarr(6)
et=dblarr(6)


for i=0, nuniq-1 do begin   ; loop over each target
  
  fileout1=outdir+'all_time_info.txt'
  fileout2=outdir+'all_temps.ps'
  fileout3=outdir+'all_medimg2.ps'
  
  filesin=file_search(indir+'*'+uname[i]+'*.sav', count=nf)
  
  for j=0, nf-1 do begin
    
    restore, filesin[j] ;roi_name, exp_num, ra_dec, jd, data1, roi_dim, xc, yc, ccd_temp, exp_time, exp_ttl, $
                          ;medcols, medimg, ndead, rightcol, simbad_radec, vmag, bmag, parlax, otype, sptype
                          
    xx=where(exp_ttl eq 10, nxx)
    print, filesin[j], nxx
    
    n=n_elements(jd)
    jd=jd-jd[0]
    jd1=jd[1:n-1]
    jdiff=jd1-jd
    gap=where(jdiff gt 0.015, ngap)
    
    stop
    
    continue
    
    
    stop
    
    st[j]=jd[0]
    et[j]=jd[n_elements(jd)-1]
                          
    ; get date range from jd 
    ; get % of time observed from number of points 
    ; get cadence
    ; open lun fileout and append data
    
    nfrm=n_elements(jd)                                                       ; SAVE = total frames
      
    single=where(exp_ttl eq 1., nsin)                                         ; SAVE
    stk=where(exp_ttl eq 10., nstk)                                           ; SAVE
    
    startdate=jd[0]       ; in julian days
    enddate=jd[nfrm-1]    ; in julian days ....convert to calendar
    
    totdur=enddate-startdate                                                  ; SAVE
    
    caldat, startdate, mon1, day1, yr1
    startdate=strtrim(day1,1)+'_'+strtrim(mon1,2)+'_'+strtrim(yr1,2)          ; SAVE
    caldat, enddate, mon2, day2, yr2
    enddate=strtrim(day2,1)+'_'+strtrim(mon2,2)+'_'+strtrim(yr2,2)            ; SAVE
    
    ; get cadence of single and stacked images
    jdsin=jd[single]
    jd2=jdsin[1:nsin-1]
    jdiff=jd2-jdsin
    cadence=median(jdiff)                                                     
    cadence_si=cadence*24.*60.*60                                                ; SAVE - in secs
    
    if nstk gt 0 then begin
    jdstk=jd[stk]
    jd2=jdstk[1:nstk-1]
    jdiff=jd2-jdstk
    cadence=median(jdiff)
    cadence_st=cadence*24.*60.*60                                                ; SAVE - in secs
    
    ; get percentage of time observed......
    ; calculate expected number of points based on exposure time and 15 min orbit every 100 mins
    tot_obs_time=(nsin*cadence_si/60.)+(nstk*cadence_st/60.)  ; mins
    possible_obs_time=totdur*24.*60./100.*15.                 ; mins
    duty_cyc=tot_obs_time/possible_obs_time*100.              ; SAVE - % number of points observed in totdur
    endif else begin
      
    cadence_st=float('Nan')
    
    ; get percentage of time observed......
    ; calculate expected number of points based on exposure time and 15 min orbit every 100 mins
    tot_obs_time=(nsin*cadence_si/60.)  ; mins
    possible_obs_time=totdur*24.*60./100.*15.                 ; mins
    duty_cyc=tot_obs_time/possible_obs_time*100.              ; SAVE - % number of points observed in totdur
    
    endelse

                              
    ; temporary to handle stacked files....
    x1=(where(finite(ccd_temp) eq 1, nx1))[0]
    x2=x1+nx1-1
    
    starttemp=ccd_temp[x1]
    endtemp=ccd_temp[x2]
    
    alljds=[alljds,jd]
    allttls=[allttls,exp_ttl]
    alltemps=[alltemps,ccd_temp]
    imgmeds=[imgmeds,medimg]
    
    med=fltarr(nfrm)
    for k=0, nfrm-1 do med[k]=median(data1[*,*,k])
    
    imgmeds2=[imgmeds2,med]
    
    xx=where(exp_ttl eq 10.)
    
    plotsym, 0, /fill, 0.7
    ;plot, ccd_temp, medimg, psym=8, color=cgcolor('black'), xrange=[15,40]
    ;oplot, ccd_temp[xx], med[xx], psym=8, color=cgcolor('purple')
    
    ;yy=where(medimg[xx] gt 1500, nyy)
 
  ;  for gg=0, nyy-1 do begin
      
  ;    data2=data1[*,*,xx[yy[gg]]]
  ;    plot_image, bytscl(data2, 20, 1000)
  ;    wait, 0.2
  ;  endfor
    

  ;  openw, lun, fileout1, /get_lun, /append
  ;  if j eq 0 then printf, lun, 'start_date', 'end_date', 'totdur_days', 'num_1s', 'num_10x', $
  ;    '1s_cadence', '10s_cadence', 'dutycyc_pcent', 'starttemp', 'endtemp', $
  ;    format='(a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15)'
  ;  printf, lun, startdate, enddate, totdur, nsin, nstk, cadence_si, cadence_st, $
  ;    duty_cyc, starttemp, endtemp, $
  ;    format='(a15,x,a15,x,d15.5,x,i15,x,i15,x,d15.5,x,d15.5,x,d15.5,x,d15.5,x,d15.5)'
  ;  free_lun, lun
   

  endfor
  
  time=alljds
  temp=alltemps
  
  nfrm=n_elements(time)
  
  caldat, time, mon, day, yr
  
  dates=strtrim(day,2)+'_'+strtrim(mon,2)+'_'+strtrim(yr,2)
  
  xx=where(allttls eq 10, complement=yy)
  
  tickloc=lonarr(6)
  obmon=['10','11','12','1','2','3']
  for k=0, 5 do tickloc[k]=(where(strtrim(mon,2) eq obmon[k]))[0]
  
  col1=cgcolor(['red', 'blue', 'green', 'orchid', 'forest green', 'purple'])
    
 ; plotsym, 0, /fill, 0.2
 ; ps_on, fileout2, xsize=18, ysize=12
 ; plot, time, temp, psym=8, color=cgcolor('black'), /ynozero, xrange=[time[0]-10., time[nfrm-1]+9], $
 ;   xstyle=1, xtickv=time[tickloc], xtickname=dates[tickloc], xticks=5, xtitle='Dates', ytitle='CCD Temp', $
 ;   charsize=0.7, title='Orion'
 ; oplot, time[xx], temp[xx], color=cgcolor('blue'), psym=8
 ; for k=0, 5 do oplot, [st[k], st[k]], [0,400], color=col1[k], thick=3
 ; for k=0, 5 do oplot, [et[k], et[k]], [0,400], color=col1[k], thick=3
 ; ps_off
 ; !p.multi=0
 ; stop
; ps_on, fileout3, xsize=18, ysize=12
  plot, time[xx], imgmeds2[xx], psym=8, color=cgcolor('black'), /ynozero, xrange=[time[0]-10., time[nfrm-1]+9], $
    xstyle=1, xtickv=time[tickloc], xtickname=dates[tickloc], xticks=5, xtitle='Dates', ytitle='Image Medians', $
    charsize=0.7, title='Orion';, yrange=[0,100]
  oplot, time[xx], imgmeds2[xx], color=cgcolor('blue'), psym=8
  for k=0, 5 do oplot, [st[k], st[k]], [0,400], color=col1[k], thick=3
  for k=0, 5 do oplot, [et[k], et[k]], [0,400], color=col1[k], thick=3
;  ps_off
  stop
  
endfor


print, 'end of program'

end

