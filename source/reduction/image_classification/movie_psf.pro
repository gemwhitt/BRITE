pro movie_psf

; make movies to show the change in target PSF towards the end of the orbits
; 
; beginning orbit, middle and near the end
; 
Compile_opt idl2
!p.background=cgcolor('white')

indir='~/BRITE/data/UB/p2/ORION/sg1/sav/'

outdir='~/BRITE/data/UB/p2/ORION/sg1/movies/'

filesin=file_search(indir+'*.sav', count=nf)

if nf eq 0 then stop

for f=0, nf-1 do begin
  
  restore, filesin[f]
  
  hdname=file_basename(filesin[f], '_sg1.sav')
  
  jd1=jd-jd[0]
  nfrm=n_elements(jd1)
  jd2=jd1[1:nfrm-1]
  
  jdiff=jd2-jd1
  
  gap=where(jdiff gt 0.015, ngap)
  gap=[-1,gap,nfrm-1]
  
  iloc1=indgen(gap[1]+1)
  iloc2=indgen(gap[501]-gap[500])+gap[500]+1
  
  ni1=n_elements(iloc1)
  
  twofrm=lonarr(70, 32, ni1)
  twofrm[0:31,0:31,*]=data1[*,*,iloc1]
  twofrm[38:69,0:31,*]=data1[*,*,iloc2[0:ni1-1]]
  twofrm[32:37,0:31,*]=1000
  
  for frm=0, ni1-1 do begin
  
    plot_image, bytscl(twofrm[*,*,frm], 20, 200), $
      title='CCD Temp = '+strtrim(ccd_temp[iloc1[frm]],2)+'          '+hdname+',          CCD Temp = '+$
      strtrim(ccd_temp[iloc2[frm]],2), $
      color=cgcolor('black'), xtickname=replicate(' ', 10), charsize=0.7
   wait, 0.5
  endfor

  continue
  
  fileout=hdname+'_sg1.mp4'  
  CD, outdir 
  mpgFilename = fileout
  
  ; Set up the video player for output.
  video = IDLffVideoWrite(mpgFilename, Format='mp4')
  framerate = 5
  framedims = [800, 500]
  stream = video.AddVideoStream(framedims[0], framedims[1], framerate)
  
  count=ni1
  
  for frm=0, count-1 do begin
    
    cgPS_Open, 'movie.ps', /Quiet
    cgDisplay, 800, 500
    plot_image, bytscl(twofrm[*,*,frm], 20, 200), title=hdname+', Exp '+strtrim(frm,2), $
      color=cgcolor('black'), xtickname=replicate(' ', 10), charsize=0.7
      
    cgPS_Close, /PNG, Width=800, /Delete_PS ; Convert to PNG file.
    
    ; Read the high-resolution PNG file.
    image = Read_PNG('movie.png')
    File_Delete, 'movie.png'
    
    ; Add the high-resolution image to the video stream.
    void = video -> Put(stream, image)

  endfor
  
  video -> Cleanup

endfor


print, 'end of program'
end