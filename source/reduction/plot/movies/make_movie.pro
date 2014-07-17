pro make_movie

; program to make movie plots to send to Rainer
; 
Compile_opt idl2

;nstk=1

target='36861'
  
indir='~/BRITE/data/UB/p1/ORION/'       ;+strtrim(nstk,2)+'stk/'

outdir='~/BRITE/data/UB/movies/ORION/'          ;+strtrim(nstk,2)+'stk/'
  
filesin=file_search(indir+'*'+target+'*.sav', count=nf)

for i=0, nf-1 do begin
                        
  obj=obj_new('IDL_Savefile', filesin[i])
  obj->restore, 'data1'
  obj->restore, 'ccd_temp'
  obj->restore, 'medimg0'
  obj->restore, 'medimg1'
  obj->restore, 'roi_name'
  obj->restore, 'xc'
  obj->restore, 'yc'
  
;  jd1=jd-jd[0]
       
  ; get x and y dimensions of data1
  xdim=(size(data1, /dim))[0]
  ydim=(size(data1, /dim))[1]
    
  ; find number of orbits - only stack within an orbit
  nfrm=n_elements(jd)
  
  fname=roi_name+'_1.mp4'
  
  CD, '~/BRITE/data/UB/movies/ORION/'
 
  mpgFilename = fname
  
  ; Set up the video player for output.
  video = IDLffVideoWrite(mpgFilename, Format='mp4')
  framerate = 10
  framedims = [600, 580]
  stream = video.AddVideoStream(framedims[0], framedims[1], framerate) 
    
  ; We are going to create high-resolution PNG images, which we will add
  ; to the video stream. We create the high-resolution PNG image from PostScript
  ; intermediate files.
  ;position = [0.1, 0.1, 0.9, 0.8]
  ;cgProgressBar = Obj_New("CGPROGRESSBAR", /Cancel, Title='Creating MPG-4 Movie...')
  ;cgProgressBar -> Start
    
  count=200
    
    FOR j=100,count-1 DO BEGIN
    
      ; Create the PostScript file.
      cgPS_Open, 'movie.ps', /Quiet
      cgDisplay, 600, 580
      cgimage, bytscl(data1[*,*,j],20, 500), title=roi_name, charsize=0.7, color=cgcolor('black')
      
      cgPS_Close, /PNG, Width=600, /Delete_PS ; Convert to PNG file.
      
      ; Read the high-resolution PNG file.
      image = Read_PNG('movie.png')
      File_Delete, 'movie.png'
      
      ; Add the high-resolution image to the video stream.
      void = video -> Put(stream, image)
      
      ; Did the user cancel the operaton?
      ;IF cgProgressBar -> CheckCancel() THEN BEGIN
      ;  ok = Dialog_Message('The user cancelled operation.')
      ;  video -> Cleanup
      ;  RETURN
      ;ENDIF ELSE cgProgressBar -> Update, (j+1)
      
    ENDFOR
    
    ; Clean up by closing the video file and writing its location.
    ;cgProgressBar -> Destroy
    video -> Cleanup
    CD, Current=currentDir
    CD, '~/BRITE/data/UB/movies/ORION'
    
    Print, 'File "' + mpgFilename + '" written to ' + currentDir + '.'
        
  endfor
 
  print, 'end of program'
end