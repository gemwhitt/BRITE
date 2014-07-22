;+
; NAME:
; FindPeaks
; PURPOSE:
; Find local maxima in 2D array
; CATEGORY:
; smei/gen/idl/toolbox
; CALLING SEQUENCE:
FUNCTION FindPeaks, Box , $
  mask  = Mask    , $
  npeak = npeak   , $
  flat  = flat    , $
  relative= relative  , $
  map   = Map   , $
  fraction= fraction  , $
  untested= untested  , $
  mindist = mindist , $
  count = count
  ; INPUTS:
  ; Box     2D array, any type
  ; OPTIONAL INPUT PARAMETERS:
  ; mask=Mask 1D or 2D array, type integer, default: no mask (includes all elements)
  ;         Elements outside the mask are ignored
  ;         1D array: list of indices defining the area in Box to be searched
  ;         for maxima.
  ;           !!! A 1-element mask MUST be entered as an array of 1-element.
  ;             A scalar is ignored, so put brackets around the input
  ;             argument if necessary)
  ;         2D array: array of same dimensions as Box with value 1 defining
  ;           the mask and value 0 outside the mask.
  ; fraction=fraction
  ;       scalar; same type as input array; default: 1.0
  ;         fraction between 0 and 1. Only this fraction of the image inside
  ;         the mask is searched for local maxima.
  ; mindist=mindist
  ;       scalar; type: any; default: 0
  ;         minimum distance between local maxima. Local maxima closer then
  ;         mindist to a higher local maximum are included with the higher
  ;         maximum.
  ; npeak=npeak scalar; type: integer
  ;         Maximum number of local maxima to be detected. By default the
  ;         whole array is searched for all maxima. If npeak is set then the
  ;         search is interrupted when npeak local maxima have been found.
  ;         In this case the 'map' array will be partially filled with the
  ;         value returned in 'untested' to indicate that this part of Box
  ;         has not been searched for maxima.
  ; flat=flat scalar; same type as input array.
  ;         In determing whether neigbour pixels are above/at/below the
  ;         central pixel of a 3x3 group of pixels, above/at/below are defined
  ;         as follows:
  ;           below:  Vneighbour < Vcenter-flat
  ;           at   :  Vcenter-flat <= Vneigbour <= Vcenter+flat
  ;           above:  Vneighbour > Vcenter+flat
  ;         In a 'noisy' map setting flat to a positive value, suppresses the
  ;         detection of lots of unwanted local maxima. (Probably the same can be
  ;         accomplished by applying the IDL smooth function to the input
  ;         argument Box before calling this function.)
  ; OUTPUTS:
  ; Peaks   1D array, type integer
  ;         Indices of local maxima inside area defined by mask
  ; OPTIONAL OUTPUT PARAMETERS:
  ; map=Map   2D array, type long integer, same dimension as Box
  ;         Provisionally divides Box into areas associated with each
  ;         local maxima. The index from the Peaks array is used
  ;         in Map to indicate the associated maximum
  ; untested=untested
  ;       value used to indicate untested areas in Box (see npeak keyword)
  ; INCLUDE:
  @compile_opt.pro    ; On error, return to caller
  ; CALLS:
  ; InitVar, IsType, boost, where_common
  ; PROCEDURE:
  ; The entire input image is searched for local maxima. This very quickly
  ; takes and intolarably long time. Usually keywords fraction and/or npeak need
  ; to be used to get results reasonably quick.
  ; MODIFICATION HISTORY:
  ; OCT-1998, Paul Hick (UCSD/CASS)
  ; JAN-2003, Paul Hick (UCSD/CASS)
  ;   Fixed problems with local maxima of equal height. Added keywords npeak,
  ;   untested and flat. Significant speedup by processing groups of up to
  ;   nine pixels (3x3 group) instead of only single pixels.
  ; FEB-2008, Paul Hick (UCSD/CASS)
  ;   Modified to deal with NaN in Box (always excluded from Mask)
  ;   Bugfix.
  ;-
  
  W = size(Box, /dim)               ; Dimensions of 2D array
  
  InitVar, flat  , 0
  InitVar, relative, /key
  InitVar, npeak   , W[0]*W[1]
  InitVar, fraction, 1
  InitVar, mindist , 0
  
  mindist2 = mindist*mindist
  
  X = lindgen(W[0])# replicate(1L,W[1])     ; X-coordinates
  Y = replicate(1L,W[0])# lindgen(W[1])     ; Y-coordinates
  
  UnTested  = -2L
  OutsideMask = -1L
  
  ; Only pixels inside mask are tested
  ; Always exclude NaNs from the mask.
  
  CASE 1 OF
    IsType(Mask, /undefined): pp = finite(Box)
    size(Mask,/n_dim) EQ 2  : pp = finite(Box) AND Mask ; Boolean mask
    ELSE: BEGIN                   ; Index array
      pp = bytarr(W)
      i  = where(finite(Box[Mask]))
      IF i[0] NE -1 THEN pp[Mask[i]] = 1
    END
  ENDCASE
  
  Map = replicate(OutsideMask,W)
  pp = where(pp)
  IF pp[0] NE -1 THEN Map[pp] = UnTested
  
  nTest = round(fraction*total(Map EQ UnTested))
  
  nTested = 0
  
  x8 = [-1, 0, 1,-1,1,-1,0,1]
  y8 = [-1,-1,-1, 0,0, 1,1,1]
  
  Pix = where(Map EQ UnTested)          ; Indices of remaining pixels to be checked
  count = 0L
  
  WHILE Pix[0] NE -1 AND n_elements(ppmax) LT npeak AND nTested LE nTest DO BEGIN
  
    HiVal = max(Box[Pix],pp,/nan)       ; Highest remaining untested value
    ; (not necessarily a maximum)
    IF relative THEN FlatVal = flat*(HiVal > 1) ELSE FlatVal = flat
    pp = Pix[pp]                ; 1D index of high value
    xp = X[pp]
    yp = Y[pp]
    
    xn = xp+x8
    yn = yp+y8
    i  = where(0 LE xn AND xn LE W[0]-1 AND 0 LE yn AND yn LE W[1]-1)
    xn = xp+x8[i]               ; X,Y indices of neighbours of high value
    yn = yp+y8[i]               ; (there will be 8 neighbours or less)
    
    i = where(Map[xn,yn] NE OutsideMask)    ; Find neighbours located inside mask
    
    IF i[0] EQ -1 THEN BEGIN          ; No neighbours inside Mask
      ; (counts as local maximum)
      ppnew = pp
      UnTestedN = -1
      
    ENDIF ELSE BEGIN              ; There are neigbours inside Mask
    
      xn = xn[i]                ; X,Y indices of neighbours inside mask
      yn = yn[i]
      
      ; Identify tested and untested neighbours
      
      UnTestedN = where(Map[xn,yn] EQ UnTested, complement=TestedN)
      
      ; Indentify neighbours below/at/above pp
      
      bb = Box[xn,yn]
      
      ; pp is always the untested pixel with the highest value
      ; (though there may be others with the same value).
      
      ; First decide whether pp is a new local maximum or not.
      
      ; If one of the neigbours lies above HiVal (Above[0] ne -1) then pp is
      ; not a local maximum. In addition, since pp is the highest untested pixel,
      ; this neighbour was tested already
      
      Above = where(bb GT HiVal+FlatVal)    ; Find neighbours above
      NewPeak = Above[0] EQ -1
      
      IF NewPeak THEN BEGIN         ; No neighbours above HiVal
      
        At = where(HiVal-FlatVal LE bb AND bb LE HiVal+FlatVal)
        NewPeak = At[0] EQ -1
        
        ; If NewPeak = 1:
        ;   all neighbours (tested and untested) lie below HiVal:
        ;   pp is a new local maximum
        
        ; If NewPeak = 0:
        ;   Some neighbours are at HiVal
        ;     All neighbours (tested and untested) lie at/below HiVal.
        ;   If none of the neigbours at HiVal have been tested already, then this
        ;   is a new plateau (several neighbouring pixels at HiVal).
        
        IF NOT NewPeak THEN NewPeak = (where_common(TestedN, At))[0] EQ -1
        
      ENDIF
      
      CASE NewPeak OF
      
        ; There are two types, which are not new peaks:
        ; 1. there is a neighbouring, already tested, pixel above HiVal
        ; 2. there is a neighbouring, already tested, pixel at HiVal
        ; Mark pp by assigning it to the nearest local maximum.
        
        0: BEGIN
        
          AtAbove = where(bb GE HiVal-FlatVal)
          TestedN = TestedN[where_common(TestedN,AtAbove)]; Tested neighbours at/above HiVal
          
          i = Map[xn[TestedN],yn[TestedN]]  ; Pick nearest local maximum
          F = min((xp-X[i])^2+(yp-Y[i])^2,i)
          i = TestedN[i]
          
          ppnew = Map[xn[i],yn[i]]
          
        END
        
        ; There are two types of new peaks:
        ; 1. all neighbours are below HiVal
        ; 2. all neighbours are at/below HiVal; all neighbours at HiVal are untested
        
        1: ppnew = pp
        
      ENDCASE
      
    ENDELSE
    
    count_old = n_elements(ppmax)
    
    CASE 1 OF
    
      ; Not a new local maximum
      
      ppnew NE pp :             ; Nothing to do
      
      ; New local maximum (ppnew = pp)
      
      mindist EQ 0: boost, ppmax, pp      ; No limits on distance between maxima
      
      IsType(ppmax,/undefined): ppmax = pp  ; Accept first local maximum unconditionally
      
      ELSE: BEGIN               ; Accumulate locations of local maxima
      
        dist2 = (xp-X[ppmax])^2+(yp-Y[ppmax])^2
        F = min(dist2,i)
        dist2 = dist2[i]          ; Distance^2 to nearest maximum
        
        ; If the new peak is closer than mindist to the nearest higher peak
        ; then cancel the maximum.
        
        CASE dist2 LT mindist2 OF
          0: boost, ppmax, pp         ; Accumulate locations of local maxima
          1: ppnew = ppmax[i]         ; Cancel new maximum
        ENDCASE
        
      END
      
    ENDCASE
    
    Map[pp] = ppnew             ; Mark position of maximum
    nTested += 1
    
    IF UnTestedN[0] NE -1 THEN BEGIN
      Map[xn[UnTestedN],yn[UnTestedN]] = ppnew
      nTested += n_elements(UnTestedN)
    ENDIF
    
    count = n_elements(ppmax)
    
    ;IF count GT count_old THEN BEGIN
    IF mindist NE 0 THEN BEGIN
      dist2 = (xp-X)^2+(yp-Y)^2
      i = where(dist2 LE mindist2 AND Map EQ UnTested)
      IF i[0] NE -1 THEN BEGIN
        Map[i] = ppnew
        nTested += n_elements(i)
      ENDIF
    ENDIF
    ;ENDIF
    
    Pix = where(Map EQ UnTested)        ; Indices of remaining pixels to be checked
    
  ENDWHILE
  
  RETURN, ppmax  &  END