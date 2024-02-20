begin
;************************************************
; read in data
;************************************************
  
  in = addfile("/gpfs/home/cd21q/cloud_top_hgt/WRF-Chem_CTH/29May2012/woLDA/2012-05-29_23:30.nc","r")
  ;in = addfile("/gpfs/home/cd21q/cloud_top_hgt/GOES_CTH/11June/CTH_DA_d02/CTH_regridded_GridSat-CONUS.goes13.2012.06.11.2200.v01.nc","r")
  ;in = addfile("/gpfs/home/cd21q/cloud_top_hgt/GOES_CTH/29May/CTH_DA_d02/CTH_regridded_GridSat-CONUS.goes13.2012.05.29.2330.v01.nc","r")


  lat1=in->lati
  lon1=in->longi
  lat1@units = "degrees_north"
  lon1@units = "degrees_east"

  a=in->CTH ;GOES13_CTH ;CTH
  Ttotal = a

  printVarSummary(Ttotal)

  max_f = max(Ttotal)
  print(max_f)

;***********************************************
; common resources
;************************************************
  wks  = gsn_open_wks("pdf","CTHs_Sprcll_2330")         ; send graphics to pdf file
  gsn_define_colormap(wks,"precip3_16lev");"GreenMagenta16")


  res                       = True         ; plot mods desired
  res@gsnDraw               = False   ; Don't draw plot (will do later)
  res@gsnFrame              = False   ; Don't advance framce  (will do later)
  res@cnFillOn              = True         ; turn on color
  res@cnFillOpacityF        = 1            ;opacity
  res@cnFillMode            = "AreaFill";"CellFill"
  res@cnLinesOn             = False        ; turn off contour lines
 
     res@tfDoNDCOverlay        = True

  res@lbLabelStride         = 1            ; skip every other label
  res@lbBoxLinesOn          = False        ; turn off box between colors

  res@gsnSpreadColors       = True         ; use full colormap
  res@gsnMaximize           = True            ; enlarge plot

  res@gsnLeftString        = "WRF-Chem 12 km ~C~ Cloud Top Heigts (km)";"WRF-Chem 12 km ~C~4-hour Pr(mm)"
  res@gsnCenterString	   = "";"with LDA";"with CTHDA~C~Cansu@FSU"
  res@gsnRightString	   = "T:2012-05-29_23:30:00";"   Init:2012-05-29_20:00:00~C~Valid:2012-05-30_00:00:00"
  res@gsnStringFontHeightF = 0.0115


  res@mpDataBaseVersion     = "MediumRes"
  res@mpDataSetName         = "Earth..4"
  res@mpProjection          = "Mercator"     ; choose map projection
; res@mpCenterLonF          =  40.           ; choose center lon
; res@mpCenterLatF          =  36.           ; choose center lat
  res@mpLimitMode           = "LatLon"       ; required
  res@mpMinLatF             =  min(lat1);34     ; min lat
  res@mpMaxLatF             =  max(lat1);43     ; max lat
  res@mpMinLonF             =  min(lon1);25     ; min lon
  res@mpMaxLonF             =  max(lon1);45.4     ; max lon

  res@mpOutlineBoundarySets     = "GeophysicalAndUSStates"; turn on states

  res@pmTickMarkDisplayMode = "Always"        ; turn on automatic tickmarks
  res@mpGridAndLimbOn       = True            ; turn on lat/lon lines
  res@mpGridMaskMode        = "MaskLand"      ; Mask grid over land.
  res@tiMainString          = "";"Annual Mean Total Precipitation"
  res@gsnAddCyclic          = False           ;
; res@trGridType            = "TriangularMesh"
;  res@sfXArray              = lon1
;  res@sfYArray              = lat1


  res@cnLevelSelectionMode     = "ExplicitLevels"

  res@cnLevels                    = (/2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14/)


;***********************************************************
; Draw maps
;***********************************************************


  txres                       = True
  txres@txPerimOn             = True
  txres@txFontHeightF         = 0.02
  txres@txBackgroundFillColor = "White";"LightGrey"
  txres@txFont                = "helvetica-bold"

  amres                       = True
  amres@amParallelPosF        = -0.41    ; This is the left edge of the plot.
  amres@amOrthogonalPosF      = 0.45     ; This is the bottom edge of the plot.
  amres@amJust                = "BottomLeft"

  lnres                       = True           ; Set some line options
  lnres@gsLineThicknessF      = 1.0            ; default is 1.0
  lnres@gsLineColor           = "Black"        ; default is black


  plot = gsn_csm_contour_map(wks,Ttotal/1000,res)    

  txid = gsn_create_text(wks, "non-LDA", txres)
  annoid = gsn_add_annotation(plot, txid, amres)
  draw(plot)
  frame(wks)

end
