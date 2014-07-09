pro read_rasters_pointing

Compile_opt idl2

devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
imagelib    ; adds system variable !IMAGE needed for plotting
  
; Program to read raster files in pointing directory, display image and observation info....
;
sat='UB'
obsdate='20130717'
indir='~/BRITE/data/'+sat+'/rasters/fits_data/'+obsdate+'*'

fitsfiles=file_search(indir+'/*.fits', count=nfits)
  
print, 'Total number of fits files in '+indir+' is:'
print, nfits

for i=0, nfits-1 do begin
  
  fitsname=file_basename(fitsfiles[i], '.fits')
  print, fitsname
   
  fits_info, fitsfiles[i], n_ext=next, extname=extname, /silent ;textout='~/Desktop/fitsinfo_'+fitsname+'.txt'

  ;print, 'Number of extensions: ',next
  ;print, 'Ext names: '+extname
  
  ; read header info
  data=mrdfits(fitsfiles[i], 0, header)
  
  ; get useful header info
  telescop=sxpar(header,'TELESCOP')
  raj2000=sxpar(header,'RAJ2000')   ; IN RADIANS - CONVERT TO DEGREES OR HHMMSS (times by 180/!pi)
  decj2000=sxpar(header,'DECJ2000') 
  exp_time=sxpar(header,'EXP_TIME')
  exp_ttl=sxpar(header,'EXP_TTL')
  exp_num=sxpar(header,'EXP_NUM')
  num_roi=sxpar(header,'NUM_ROI')
  jd_obs=sxpar(header,'JD-OBS')
  date_obs=sxpar(header,'DATE-OBS')
  hjd=sxpar(header,'HJD')
  hjd_corr=sxpar(header,'HJD-CORR')
  
  ; print some header info
  print, 'Exposure number is: ', exp_num
  ;print, 'Julian date is: ',jd_obs
  print, 'Date OBS is :', date_obs
  print, 'RA (degrees) is ',raj2000*(180./!pi)
  print, 'DEC (degrees) is ',decj2000*(180./!pi)
  ; open fits file and read extensions starting at 1
  lun=fxposit(fitsfiles[i], 1)
  
  extension=1
  
  repeat begin
  
    data=mrdfits(lun, 0, hdr, status=status)
    
    ; get header info for this ROI
    naxis1=sxpar(hdr, 'NAXIS1')
    naxis2=sxpar(hdr, 'NAXIS2')
    extname=sxpar(hdr, 'EXTNAME')
    roi_indx=sxpar(hdr, 'ROI_INDX')
    x1=sxpar(hdr, 'ROI_X1')
    x2=sxpar(hdr, 'ROI_X2')
    y1=sxpar(hdr, 'ROI_Y1')
    y2=sxpar(hdr, 'ROI_Y2')
    ra=sxpar(hdr, 'NOM_RA')
    dec=sxpar(hdr, 'NOM_DEC')
    
    print, extname
    print, 'Ra=',ra, ' Dec=',dec
    
    ; calculate center of FOV
    xc=((x2-x1)/2.)+x1
    yc=((y2-y1)/2.)+y1
    
    print, 'center of ROI is ', xc, yc
    
    window, 0, xsize=500, ysize=500, xpos=1500, ypos=500
    plot_image, data  ;, keep_aspect_ratio=1
    
    ; look for peaks in the data
    mthr=100
    spr=5
    pks = brite_findpeaks(data,minthresh=mthr,spread=spr)
    
    ;find brightest peak closest to image center
    ; first compute distance between peaks and image center (35,35)
    dist1=sqrt((pks[0,*]-35.)^2 + (pks[1,*]-35)^2)
    ; then sort in increasing order
    order=sort(dist1)
    ; select first one
    approx_loc=fltarr(2,1)
    approx_loc[0,0]=pks[0,order[0]] ; approx x location
    approx_loc[1,0]=pks[1,order[0]] ; approx y location
    
   ; print, 'pks: ', size(pks, /dim)
    print, 'approx center of target is ', approx_loc
    
    ; replot image centered on target'
    data_target=data[(approx_loc[0,0])[0]-10:(approx_loc[0,0])[0]+10,(approx_loc[1,0])[0]-10:(approx_loc[1,0])[0]+10]
    plot_image, data_target  ;, keep_aspect_ratio=1
    
    ; try different centroid functions to get center of peak of central star
    
    ; centroid function
    centroid, data, xcen1, ycen1, xy_peak=[(approx_loc[0,0])[0],(approx_loc[1,0])[0]], fwhm=2;, /peak_loc
    centroid, data, xcen2, ycen2, xy_peak=[35,35], fwhm=2;, /peak_loc
    
    print, 'result of centroid1 is: xcen1=',xcen1,' ycen1=',ycen1
    
    print, 'result of centroid2 is: xcen2=',xcen2,' ycen2=',ycen2
    
    imstars=steve_centroid2(data,approx_loc,width=5,fwhm=2)
      
    print, 'result of steve_centroid is: ',imstars
    
    ; plot window centered on target with results of each centroid overplot
    window, 2, xsize=500, ysize=500, xpos=2500, ypos=-500
    plot_image,data[approx_loc[0,*]-10:approx_loc[0,*]+10,approx_loc[1,*]-10:approx_loc[1,*]+10], origin=[approx_loc[0,*]-10,approx_loc[1,*]-10]
    oplot, [xcen1], [ycen1], color=cgcolor('red'), psym=4, symsize=2  
    oplot, [xcen2], [ycen2], color=cgcolor('green'), psym=5, symsize=2  
    oplot, [imstars[0,0]], [imstars[1,0]], color=cgcolor('blue'), psym=2, symsize=2  
    
    stop
 
    extension=extension+1
  endrep until extension eq num_roi+1
  
  free_lun, lun
  stop
  
endfor

stop
end
