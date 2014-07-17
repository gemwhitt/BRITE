pro measure_systrend

  ; Modified from calculate_trend_orbit.pro 
  
  Compile_opt idl2
  !p.background=cgcolor('white')
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; input directory
  indir='~/BRITE/data/UB/p3/'
  
  filesin=file_search(indir+'*_p3.sav', count=nsav)
  
  obj=obj_new('IDL_Savefile', filesin[0])
  obj->restore, 'jd'
  obj->restore, 'ccd_temp'
  
  ; define new variables for saving in output file - normalized time and header-board temps
  jd1=jd-jd[0]
  header_temps=ccd_temp
  
  ; number of images/frames in this set of files
  nfrm=n_elements(jd)
  
  ; define new variables for saving in output file 
  roi_out=fltarr(nsav,nfrm) ; the mean value of pixels outside the target radius
  roi_loc=fltarr(2,nsav)    ; location of each roi on the CCD
  hdname=strarr(nsav)
  
  for bb=0, nsav-1 do begin
  
    obj=obj_new('IDL_Savefile', filesin[bb])
    obj->restore, 'data1'
    obj->restore, 'jd'
    obj->restore, 'xc'
    obj->restore, 'yc'
    obj->restore, 'roi_name'
    
    hdname[bb]=roi_name[0]
    
    roi_loc[0,bb]=xc
    roi_loc[1,bb]=yc
    
    for rr=0, nfrm-1 do begin
    
      dat=reform(data1[*,*,rr])
      rbin=4
      dat2=rebin_data(dat, rbin)

      if rr eq 0 then begin ; determine index locations of outside pixels
      
      ; get roi x and y center coords in binned pixels
      roi_xc=float((size(dat2, /dim))[0])/2.
      roi_yc=float((size(dat2, /dim))[1])/2.
      
      ; find approx center of PSF
      mthr=1.5*robust_mean(dat2,3)
      spr=2*rbin
      pks=brite_findpeaks(dat2, minthresh=mthr, spread=spr)
      
      if pks[0] eq -999.9 then goto, next_frame
      
      npks=n_elements(pks)/3.
      
      if npks gt 1 then pks=pks[*,where(pks[2,*] eq max(pks[2,*]))]
      
      pks1=round(pks)
      
      ; determine which quadrent target is in on the CCD or if it is in the center
      xd=pks[0]-roi_xc
      yd=pks[1]-roi_yc
      
      if abs(xd) lt 4.*rbin AND abs(yd) lt 4.*rbin then begin ; target is roughly in center
          ; define a region around the perimeter
          
          temp_dat=dat2
          temp_dat[roi_xc-(12.*rbin):roi_xc+(12.*rbin),roi_yc-(12.*rbin):roi_yc+(12.*rbin)]=-9999
          
          xx=where(temp_dat ne -9999, nxx)
          
          x_loc=(array_indices(dat2, xx))[0,*]
          y_loc=(array_indices(dat2, xx))[1,*]
          
          plot_image, dat2
          oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
          
         wait, 1
          
      endif else begin  ; target is within a quadrent
          ; define a region based in the other 3 quads
          
          temp_dat=dat2
          
          ; define 4 quads - lower-left, right and upper-left, right
          if xd lt 0. and yd lt 0. then temp_dat[0:roi_xc+(7.*rbin),0:roi_yc+(7.*rbin)]=-9999 ; lower-left
          if xd gt 0. and yd gt 0. then temp_dat[roi_xc-(7.*rbin):(roi_xc*2)-1, 0:roi_yc+(7.*rbin)]=-9999  ; lower-right
          if xd lt 0. and yd gt 0. then temp_dat[0:roi_xc+(7.*rbin),roi_yc-(7*rbin):(roi_yc*2)-1]=-9999  ; upper-left
          if xd gt 0. and yd gt 0. then temp_dat[roi_xc-(7.*rbin):(roi_xc*2)-1,roi_yc-(7*rbin):(roi_yc*2)-1]=-9999  ; upper-right
          
          xx=where(temp_dat ne -9999, nxx)
          
          x_loc=(array_indices(dat2, xx))[0,*]
          y_loc=(array_indices(dat2, xx))[1,*]
          
          plot_image, dat2
          oplot, x_loc, y_loc, color=cgcolor('orange'), psym=2
            
wait, 1

      endelse
      
      endif
      
      ; calculate the median in the in- and out- of target pixels
      roi_out[bb,rr]=median(dat2[xx])
      ;stop

      next_frame:
    endfor  ; end loop over individual observations
    
  endfor
  
  ; save roi_out with jd1
  outfile='/Users/gemmawhittaker/BRITE/reduction/orbit_trends/roi_med_all/roi_med_all_p3.sav'
  save, filename=outfile, roi_out, jd1, header_temps, roi_loc, hdname
  
  print, 'End of Program'
  
end