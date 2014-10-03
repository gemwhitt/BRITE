pro master_program2

  ; Program to run all programs successively
  ;
  Compile_opt idl2
  
  sat='TOR'

field='CENTAURUS'

target=['HD127973','HD129056']
target='*'
  
  print, 'running aper_phot_sav'
  
  aper_phot_sav, sat, field, target ; p5 to aper_lc_sav
  ;sos_phott, sat, field, target ; p5 to aper_lc_sav
  
  print, 'running best_aperture'
  
  best_aperture, sat, field, target ; aper_lc_sav to aper_lc_best
  
  print, 'running int_pix_corr'
  
  int_pix_corr, sat, field, target  ; aper_lc_best to aper_lc_intpix
  
  print, 'running final_sigclip2'
  
  final_sigclip2, sat, field, target  ; aper_lc_intpix to aper_lc_final
  
  print, 'concatenating files'
  
  concat_lcs, sat, field, target ; aper_lc_final updated
  
  ;make_lc_plot, sat, field
  
  ;make_stats_plot, sat, field
  
  ;
  ;
  print, 'end of program'
  print, 'analyse light curves'
end