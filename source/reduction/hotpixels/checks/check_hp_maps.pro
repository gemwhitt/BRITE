pro check_hp_maps

; program to check results of HP removal - i.e. view images from p2 files

Compile_opt idl2
;; FOR PLOTTING
devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
!p.background=cgcolor('white')

indir='~/BRITE/data/UB/p2/'

filesin=file_search(indir+'Orion*.sav', count=nsav)

for i=45, nsav-1 do begin
  
  print, file_basename(filesin[i], '_p2.sav')
  
  restore, filesin[i] ;roi_name, exp_num, ra_dec1, jd, data1, roi_dim, xc, yc, ccd_temp, simbad_radec, $
                      ;simbad_mags, parlax, otype, sptype, exp_time, exp_ttl
                      
  stacked=exp_ttl/exp_time
  print, stacked
  
  data1=data1/stacked
  
  print, exp_time
  
  nfrm=(size(data1, /dim))[2]
  
  ; choose to look at 2% of all images
  for img=0, nfrm-1 do if img mod 20. eq 0 then begin
    
    plot_image, bytscl(data1[*,*,img], 50, 500)
  
  wait, 1 
  endif
  
endfor


stop
print, 'end of program'
end