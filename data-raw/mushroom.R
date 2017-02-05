example2 <- data.frame(Species = c("arorae","arvenis","benesi","bernardii","bisporus","bitorquis","califorinus","campestris","comtulus",
                                   "cupreoBrunneus","diminutives","fuseoFibrillosus","fuscovelatus","hondensis","lilaceps","micromegathus",
                                   "praeclaresquamosus","pattersonae","perobscurus","semotus","silvicola","subrutilescens","xanthodermus"),
                       PileusWidth_L = c(3,6,4,7,5,5,4,5,2.5,2.5,1.5,4,3.5,7,8,2.5,7,5,8,2,6,6,5),
                       PileusWidth_U = c(8,21,8,16,12,15,11,10,4,6,2.5,15,8,14,20,4,19,15,12,6,12,12,17),
                       StipeLength_L = c(4,4,5,4,2,4,3,3,3,1.5,3,4,4,8,9,2.5,8,6,6,3,6,6,4),
                       StipeLength_U = c(9,14,11,7,5,10,7,6,5,3.5,6,15,10,14,19,4.5,15,15,12,7,12,16,14),
                       StipeThickness_L = c(0.5,1,1,3,1.5,2,0.4,1,0.4,1,0.25,1.5,1,1.5,3,0.1,2,2.5,1.5,0.4,1.5,1,1),
                       StipeThickness_U = c(2.5,3.5,2,4.5,2.5,4.0,1.0,2.0,0.7,1.5,0.35,2.5,2,2.5,5,0.7,3.5,3.5,2,0.8,2,2,3.5))

save(example2, file = "data/example2.rdata")
