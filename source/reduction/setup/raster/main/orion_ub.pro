pro orion_ub

Compile_opt idl2
  
indir='/Users/gemmawhittaker/BRITE/data/UB/ORION/roi_raw/sav/'
    
txtout='/Users/gemmawhittaker/BRITE/results/ba_ub/ORION/ub_temps.txt'
  
filesin=file_search(indir+'HD31237*.sav', count=nf)
  
restore, filesin  ;roi_name, exp_num, ra_dec, jd, roi_loc, ccd_temp, exp_time, exp_ttl, $
                      ;simbad_radec, vmag, bmag, parlax, otype, sptype, medimg0, ndead, nsat, medcol1, data1
    
nimg=n_elements(jd)
    
ccd_temp2=fltarr(nimg)
for im=0, nimg-1 do ccd_temp2[im]=average(ccd_temp[*,im])
    
temps=ccd_temp2

; find jd nearest to 2456628.428987
jdiff=abs(jd-double(2456628.428987))

xx=min(jdiff, ii)
stop
jd0=jd[ii:nimg-1]
temps=temps[ii:nimg-1]

nimg=n_elements(jd0)
    
openw, lun, txtout, /get_lun
for ii=0, nimg-1 do printf, lun, jd0[ii], temps[ii], format='(d14.6,x,f7.3)'
free_lun, lun
   
  stop
end

