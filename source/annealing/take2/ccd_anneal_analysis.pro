pro ccd_anneal_analysis

  ; program to analyse the residuals in the columns and pixels after annealing
  
  Compile_opt idl2
  
  plotout='y'
  
  ; FOR PLOTTING
  devicelib   ; adds system variables !BCOLOR and !ASPECT needed for plotting
  imagelib    ; adds system variable !IMAGE needed for plotting
  !p.background=cgcolor('white')
  ;
  ; input directory
  indir='/Users/gemmawhittaker/BRITE/data/anneal_test_data/sav_files/'
  prefile=file_search(indir+'pre.sav')
  postfile=file_search(indir+'post.sav')
  
  ; save results to....
  outdir='~/BRITE/results/annealing/date_8nov13/'
  
  readtimes=[0.06,1.,10.]
  ntime=n_elements(readtimes)
    
    for i=1, 1 do begin ;ntime-1 do begin
    
      restore, prefile  ;fdate, temp, data, readtime  ;DATA=LONG_Array[4048, 2672, 3]
      
      pre_dat=data[*,*,i]
      pre_temp=temp[i]
     
      ; now restore post-data
      restore, postfile  ;fdate, temp, data, readtime  ;DATA=LONG_Array[4048, 2672, 3]
      
      post_dat=data[*,*,i]
      post_temp=temp[i]
      
      ; identify zero-valued pixels in pre and post data?
      pre_z=where(pre_dat eq 0, nprez)
      post_z=where(post_dat eq 0, npostz)
      
      print, 'pre ',nprez
      print, 'post ', npostz
      
      
      
     stop
 
    endfor
  
 
  stop
  print, 'End of Program'
end