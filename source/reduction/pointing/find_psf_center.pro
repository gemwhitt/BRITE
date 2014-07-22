pro find_psf_center

; First attempt - modified from read_rasters_pointing.pro

; Program to find the center of the source PSF .....
; First attempt - using brite_findpeaks to get approx location of source, then a modified version of steven centroid, which uses the centroid fn to get the PSF center
; Output arrays include:
; exposure number, julian date
; center of PSF 
; center of PSF wrt to center of image
; ra and dec 
; 
; Calls: brite_findpeaks and steve_centroid2 to find approximate location of PSF, then find accurate location
; 
Compile_opt idl2

!p.background=cgcolor('black')

results='hc' ;  if want hardcopy of results, pc to print to screen 
  
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  
  ; Program to read raster files in pointing directory, display image and observation info....
  ;
  sat='UB'
  obsdate='20130717'
  indir='~/BRITE/data/'+sat+'/rasters/fits_data/'+obsdate+'_Rasters_*'
  outdir='~/BRITE/results/'+sat+'/pointing'
  
  fitsfiles=file_search(indir+'/*.fits', count=nfits)
  
  print, 'Total number of fits files in '+indir+' is:'
  print, nfits
  
  for i=0, nfits-1 do begin
  
    fitsname=file_basename(fitsfiles[i], '.fits')
    print, fitsname
    
    fits_info, fitsfiles[i], n_ext=nexts, extname=extnames, /silent ;textout='~/Desktop/fitsinfo_'+fitsname+'.txt'
    
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
    ;print, 'Date OBS is :', date_obs
    print, 'RA (degrees) is ',raj2000*(180./!pi), ' DEC (degrees) is ',decj2000*(180./!pi)
    
    
    ; open fits file and read extensions starting at 1
    lun=fxposit(fitsfiles[i], 1)
    
    extension=1
    
    repeat begin
    
      data=mrdfits(lun, 0, hdr, status=status)
      
      ; get header info for this ROI
      naxis1=float(sxpar(hdr, 'NAXIS1'))  ; use this variable to calculate center of each raster
      naxis2=float(sxpar(hdr, 'NAXIS2'))
      extname=sxpar(hdr, 'EXTNAME')
      roi_indx=sxpar(hdr, 'ROI_INDX')
      x1=sxpar(hdr, 'ROI_X1')
      x2=sxpar(hdr, 'ROI_X2')
      y1=sxpar(hdr, 'ROI_Y1')
      y2=sxpar(hdr, 'ROI_Y2')
      ra=sxpar(hdr, 'NOM_RA')
      dec=sxpar(hdr, 'NOM_DEC')
      
      outfile=outdir+'/datfiles2/'+strtrim(extname,2)+'_'+strtrim(roi_indx,2)+'_'+obsdate+'.dat'
      
      print, extname
      print, 'Ra=',ra, ' Dec=',dec
      
      ; calculate center of FOV
      xc=((x2-x1)/2.)+x1
      yc=((y2-y1)/2.)+y1
      
      ;print, 'center of ROI is ', xc, yc
      
      if results eq 'ps' then begin
        window, 0, xsize=500, ysize=500, xpos=1500, ypos=500
        plot_image, data  ;, keep_aspect_ratio=1
      endif
      
      ; look for peaks in the data
      mthr=100
      spr=5
      pks = brite_findpeaks(data,minthresh=mthr,spread=spr)
      
      ;find brightest peak closest to image center
      ; first compute distance between peaks and image center
      dist1=sqrt((pks[0,*]-25.)^2 + (pks[1,*]-25.)^2)
      ; then sort in increasing order
      order=sort(dist1)
      ; select first one
      approx_loc=fltarr(2,1)
      approx_loc[0,0]=pks[0,order[0]] ; approx x location
      approx_loc[1,0]=pks[1,order[0]] ; approx y location
      
      ; print, 'pks: ', size(pks, /dim)
      ;print, 'approx center of target is ', approx_loc
      stop
      if results eq 'ps' then begin
        ; replot image centered on target'
        data_target=data[approx_loc[0,0]-10:approx_loc[0,0]+10,approx_loc[1,0]-10:approx_loc[1,0]+10]
        window, 1, xsize=500, ysize=500, xpos=1500, ypos=-500
        plot_image, data_target  ;, keep_aspect_ratio=1
      endif
      
      
      imstars2=steve_centroid2(data,approx_loc,width=5,fwhm=2)
      
      ;actual location of centroid needs to have x1 and y1 added onto imstars
      xcentroid=imstars2[0,0]+x1
      ycentroid=imstars2[1,0]+y1
      
      if results eq 'ps' then begin
        print, 'result of steve_centroid is: ',imstars2
        
        ; plot window centered on target with results of each centroid overplot
        window, 2, xsize=500, ysize=500, xpos=2500, ypos=-500
        plot_image, data[25-10:25+10,25-10:25+10], origin=[25-10,25-10]
        oplot, [xcentroid], [ycentroid], color=cgcolor('blue'), psym=2, symsize=2
        oplot, [xc], [yc], psym=1, color=cgcolor('red'), symsize=2
      endif
      
      ;if results eq 'hc' then begin
      ;  imagefile=outdir+'/images/'+strtrim(extname,2)+'_'+strtrim(roi_indx,2)+'_'+strtrim(exp_num,2)+'_'+obsdate+'.ps'
      ;  ps_on, imagefile, xsize=17, ysize=17
      ;  plot_image,data[xc-10:xc+10,yc-10:yc+10], origin=[xc-10,yc-10], $
      ;    title=strtrim(extname,2)+', Exp num: '+strtrim(exp_num,2), xtitle='X-Pixel', ytitle='Y-Pixel', charsize=0.6
      ;  oplot, [xcentroid], [ycentroid], color=cgcolor('blue'), psym=2, symsize=2
      ;  oplot, [xc], [yc], psym=1, color=cgcolor('red'), symsize=2
      ;  ps_off
        
        ;spawn, 'open '+imagefile+' &'
      ;endif
     
      ; calculate xdiff, ydiff, xydiff (or absolute difference)
      x_diff=xc-xcentroid
      y_diff=yc-ycentroid
      xy_diff=sqrt((xc-xcentroid)^2+(yc-ycentroid)^2)
      
      if results eq 'hc' then begin
        ; write out results
        openw, lun1, outfile, /get_lun, /append
        printf, lun1, sat, exp_num, jd_obs, ra, dec, xc, yc, xcentroid, ycentroid, x_diff, y_diff, xy_diff, $
          format='(a2, x, i4, x, d14.6, x, d14.6, x, d14.6, x, f7.2, x, f7.2, x, f7.2, x, f7.2, x, f7.3, x, f7.3, x, f7.3)'
        free_lun, lun1
      endif
     
      if results eq 'ps' then stop
            
      extension=extension+1
      endrep until extension eq num_roi+1
    
      free_lun, lun
    
  endfor
  
print, 'End of Program'
print, 'Now run combine_pdf.pro to combine files for each ROI and produce 2x6 PDF'

end
