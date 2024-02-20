begin
;************************************************
; read in data
;************************************************

;do AAAA=12,22,1
 AAAA=20
 BBBB=00;AAAA+1
  in = addfile("/gpfs/research/eoas/cd21q/WRF-Chem_outputs/29May2012/adjustedPlevs.entrain.NAM12/wLDA/rawWRF/wrfout_d02_2012-05-29_"+AAAA+":00:00","r")  ;first
  in2 = addfile("/gpfs/research/eoas/cd21q/WRF-Chem_outputs/29May2012/adjustedPlevs.entrain.NAM12/wLDA/rawWRF/wrfout_d02_2012-05-30_0"+BBBB+":00:00","r") ;last

  lat1=in->XLAT
  lon1=in->XLONG
  lat1@units = "degrees_north"
  lon1@units = "degrees_east"

  a=in->RAINNC
  b=in->RAINC
  c=in->RAINSH  

  Tfirst=a+b;a;+b;+c

  aa=in2->RAINNC
  bb=in2->RAINC
  cc=in2->RAINSH

  Tlast=aa+bb;aa;+bb;+cc


  Ttotal=Tlast(0,:,:)-Tfirst(0,:,:);dim_sum_n(Tlast(0:0,:,:),0)-dim_sum_n(Tfirst(0:0,:,:),0)

;  Ttotal@lat2d = lat1
;  Ttotal@lon2d = lon1

  Ttotal@lat=lat1
  Ttotal@lon=lon1

  printVarSummary(Tfirst)
  printVarSummary(Tlast)
  printVarSummary(Ttotal)

  max_f = max(Ttotal)
  print(max_f)


;***********************************************
; common resources
;************************************************
  wks  = gsn_open_wks("pdf","total_29may_wLDA"+BBBB)         ; send graphics to PNG file
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

  res@gsnLeftString        = "WRF-Chem 12 km ~C~1-hour Pr(mm)"
  res@gsnCenterString	   = "   with LDA~C~Cansu@FSU";"without LDA~C~Cansu@FSU";"   with LDA~C~Cansu@FSU"
  res@gsnRightString	   = "   Init:2012-06-11_"+AAAA+":00:00~C~Valid:2012-06-12_0"+BBBB+":00:00"
  res@gsnStringFontHeightF = 0.0115


  res@mpDataBaseVersion     = "MediumRes"
  res@mpDataSetName         = "Earth..4"
  res@mpProjection          = "Mercator"     ; choose map projection
; res@mpCenterLonF          =  40.           ; choose center lon
; res@mpCenterLatF          =  36.           ; choose center lat
  res@mpLimitMode           = "LatLon"       ; required
  res@mpMinLatF             =  min(lat1);lat1(0,5,30);min(lat1);34     ; min lat
  res@mpMaxLatF             =  max(lat1);lat1(0,70,110);max(lat1);43     ; max lat
  res@mpMinLonF             =  min(lon1);lon1(0,5,30);min(lon1);25     ; min lon
  res@mpMaxLonF             =  max(lon1);lon1(0,70,110);max(lon1);45.4     ; max lon

print(min(lat1))
print(max(lat1))
print(min(lon1))
print(max(lon1))

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

  res@cnLevels                    = (/0.1,5.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,70.0,80.0/)



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


  plot = gsn_csm_contour_map(wks,Ttotal,res)    

  txid = gsn_create_text(wks, "LDA", txres)
  annoid = gsn_add_annotation(plot, txid, amres)
  draw(plot)
  frame(wks)

;end do
end
