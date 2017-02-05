example4 <- data.frame(model = c("AstonMartin","AudiA6","AudiA8","BMW7","Ferrari","HondaNSR","MecedesC","Porsche"),
                       price_L = c(260.5,68.2,123.8,104.9,240.3,205.2,55.9,147.7),
                       price_U = c(460,140.3,171.4,276.8,391.7,215.2,115.2,246.4),
                       maxvelocity_L = c(298,216,232,228,295,260,210,280),
                       maxvelocity_U = c(306,250,250,240,298,270,250,305),
                       acceltime_L = c(4.7,6.7,5.4,7,4.5,5.7,5.2,4.2),
                       acceltime_U = c(5,9.7,10.1,8.6,5.2,6.5,11,5.2),
                       cyl_L = c(5935,1781,2771,2793,3586,2977,1998,3387),
                       cyl_U = c(5935,4172,4172,5397,5474,3179,3199,3600))

save(example4, file = "data/example5.rdata")
