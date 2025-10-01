# the numbers in this script come from the output of generate_upset_input.py
# the output is the essential genes intersections upset plot

#Load and Install packages if necessary
if (!require("UpSetR")) install.packages("UpSetR")

library(UpSetR)

input <- c(
  "S.Equi_sb_equi"                       = 401,
  "S.Iniae"                              = 476,
  "S.Uberis"                             = 384,
  "S.Pneumoniae"                         = 536,
  "S.Suis"                               = 342,
  "S.Agalactiae"                         = 160,
  
  # 2-way intersections
  "S.Equi_sb_equi&S.Iniae"               = 247,  # non-sig
  "S.Equi_sb_equi&S.Uberis"              = 243,
  "S.Equi_sb_equi&S.Pneumoniae"          = 226,
  "S.Equi_sb_equi&S.Suis"                = 216,
  "S.Equi_sb_equi&S.Agalactiae"          = 124,
  "S.Iniae&S.Uberis"                     = 223,  # non-sig
  "S.Iniae&S.Pneumoniae"                 = 197,  # non-sig
  "S.Iniae&S.Suis"                       = 187,  # non-sig
  "S.Iniae&S.Agalactiae"                 = 113,
  "S.Uberis&S.Pneumoniae"                = 199,
  "S.Uberis&S.Suis"                      = 172,
  "S.Uberis&S.Agalactiae"                = 109,
  "S.Pneumoniae&S.Suis"                  = 183,
  "S.Pneumoniae&S.Agalactiae"            = 100,
  "S.Suis&S.Agalactiae"                  = 93,
  
  # 3-way intersections
  "S.Equi_sb_equi&S.Iniae&S.Uberis"      = 192,
  "S.Equi_sb_equi&S.Iniae&S.Pneumoniae"  = 171,  # non-sig
  "S.Equi_sb_equi&S.Iniae&S.Suis"        = 170,  # non-sig
  "S.Equi_sb_equi&S.Iniae&S.Agalactiae"  = 107,
  "S.Equi_sb_equi&S.Uberis&S.Pneumoniae" = 179,
  "S.Equi_sb_equi&S.Uberis&S.Suis"       = 164,
  "S.Equi_sb_equi&S.Uberis&S.Agalactiae" = 106,
  "S.Equi_sb_equi&S.Pneumoniae&S.Suis"   = 158,
  "S.Equi_sb_equi&S.Pneumoniae&S.Agalactiae" = 96,
  "S.Equi_sb_equi&S.Suis&S.Agalactiae"   = 91,
  "S.Iniae&S.Uberis&S.Pneumoniae"        = 157,  # non-sig
  "S.Iniae&S.Uberis&S.Suis"              = 145,  # non-sig
  "S.Iniae&S.Uberis&S.Agalactiae"        = 96,   # non-sig
  "S.Iniae&S.Pneumoniae&S.Suis"          = 138,  # non-sig
  "S.Iniae&S.Pneumoniae&S.Agalactiae"    = 86,
  "S.Iniae&S.Suis&S.Agalactiae"          = 82,
  "S.Uberis&S.Pneumoniae&S.Suis"         = 136,
  "S.Uberis&S.Pneumoniae&S.Agalactiae"   = 90,
  "S.Uberis&S.Suis&S.Agalactiae"         = 79,
  "S.Pneumoniae&S.Suis&S.Agalactiae"     = 77,
  
  # 4-way intersections
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Pneumoniae"        = 147,
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Suis"              = 140,
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Agalactiae"        = 93,
  "S.Equi_sb_equi&S.Iniae&S.Pneumoniae&S.Suis"          = 129,
  "S.Equi_sb_equi&S.Iniae&S.Pneumoniae&S.Agalactiae"    = 83,
  "S.Equi_sb_equi&S.Iniae&S.Suis&S.Agalactiae"          = 80,
  "S.Equi_sb_equi&S.Uberis&S.Pneumoniae&S.Suis"         = 130,
  "S.Equi_sb_equi&S.Uberis&S.Pneumoniae&S.Agalactiae"   = 88,
  "S.Equi_sb_equi&S.Uberis&S.Suis&S.Agalactiae"         = 80,
  "S.Equi_sb_equi&S.Pneumoniae&S.Suis&S.Agalactiae"     = 76,
  "S.Iniae&S.Uberis&S.Pneumoniae&S.Suis"                = 116,  # non-sig
  "S.Iniae&S.Uberis&S.Pneumoniae&S.Agalactiae"          = 79,   # non-sig
  "S.Iniae&S.Uberis&S.Suis&S.Agalactiae"                = 71,   # non-sig
  "S.Iniae&S.Pneumoniae&S.Suis&S.Agalactiae"            = 68,   # non-sig
  "S.Uberis&S.Pneumoniae&S.Suis&S.Agalactiae"           = 70,   # non-sig
  
  # 5-way intersections
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Pneumoniae&S.Suis"             = 113,
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Pneumoniae&S.Agalactiae"       = 77,
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Suis&S.Agalactiae"             = 71,
  "S.Equi_sb_equi&S.Iniae&S.Pneumoniae&S.Suis&S.Agalactiae"         = 67,
  "S.Equi_sb_equi&S.Uberis&S.Pneumoniae&S.Suis&S.Agalactiae"        = 70,
  "S.Iniae&S.Uberis&S.Pneumoniae&S.Suis&S.Agalactiae"               = 62,
  
  # 6-way core
  "S.Equi_sb_equi&S.Iniae&S.Uberis&S.Pneumoniae&S.Suis&S.Agalactiae" = 62
)


# Plot
upset(fromExpression(input), 
      nintersects = 70, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = c(1.1, 1.1, 1.1, 1.1, 1.1, 0.8), 
      point.size = 2.8, 
      line.size = 1,
      set_size.show = TRUE,
      mainbar.y.label = "Essentialome", sets.x.label = "total number of gene intersections",
)

#make bar plot
