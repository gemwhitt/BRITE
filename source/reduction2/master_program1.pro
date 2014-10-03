pro master_program1

; Program to run all programs successively
; 
Compile_opt idl2

sat='TOR'

field='CENTAURUS'

target=['HD127973','HD129056']
target='*'

print, 'running hps1b'
;;
;hps1b, sat, field, target  ; p2 to p3 files
;
print, 'running hps2'

;hps2, sat, field, target  ; p3 to p4 files

print, 'running get_psf_boundary'

;get_psf_boundary, sat, field, target  ; p4 updated

print, 'running get_com'

;get_com, sat, field, target ; p4 updated

print, 'running psf_aper'

psf_aper, sat, field, target  ; p4 to p5

;
;
print, 'end of master_program1'
print, 'run photometry programs - master_program2'
end