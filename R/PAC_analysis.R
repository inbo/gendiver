
pac26005.filters = function(add_extra=F){
  PAC_list = list()

  # Asia
  PAC_list$E2025STF198 = c(
    "Tanichthys albonubes",
    "Helostoma temminkii",
    "Trichopodus leerii",
    "Microdevario kubotai",
    "Pethia padamya")

  # ZM 1
  PAC_list$E2025STF185 = c(
    "Hyphessobrycon amandae",
    "Nannostomus beckfordi",
    "Gymnocorymbus ternetzi",
    "Jordanella floridae")

  # ZM2
  PAC_list$E2025STF1368 = c(
    "Hoplisoma sterbai",
    "Osteogaster rabauti")

  # Tanganyika
  PAC_list$E2025STF200 = c(
    "Julidochromis_unclassified",
    "Neolamprologus brevis")

  # Malawi
  PAC_list$E2025STF196 = c(
    "Aulonocara stuartgranti",
    "Cichlidae_unclassified")

  if (add_extra){
    PAC_list$extra = c("Xiphophorus hellerii")
  }

  return(PAC_list)
}


pac26005.colors = function(add_extra=F){
  # fixed colors
  mypals = list()

  # Add a priori
  # ZM1
  mypals$species["Hyphessobrycon amandae"] = "lightgreen"
  mypals$species["Nannostomus beckfordi"] = "green"
  mypals$species["Gymnocorymbus ternetzi"] = "green4"
  mypals$species["Jordanella floridae"] = "green3"

  # ZM2
  mypals$species["Hoplisoma sterbai"] = "yellow"
  mypals$species["Osteogaster rabauti"] = "yellow3"

  # Malawi
  mypals$species["Aulonocara stuartgranti"] = "lightblue1"
  mypals$species["Cichlidae_unclassified"] = "lightblue3"

  # Tanganyika
  mypals$species["Julidochromis_unclassified"] = "blue2"
  mypals$species["Neolamprologus brevis"] = "blue4"

  # Azie
  mypals$species["Helostoma temminkii"] = "salmon3"
  mypals$species["Microdevario kubotai"] = "red3"
  mypals$species["Pethia padamya"] = "salmon"
  mypals$species["Tanichthys albonubes"] = "red"
  mypals$species["Trichopodus leerii"] = "red4"


  # EXTRA
  if (add_extra){
    mypals$species["Xiphophorus hellerii"] = "purple"
  }


  return(mypals$species)
}


pac26005.palette = function(x, add_extra=F){
  # init list
  pacpal = list()
  pacpal[x] = NA
  # add PAC harcoded colors
  palcols = pac26005.colors(add_extra)
  pacpal[names(palcols)] = palcols
  # add greys for non-PAC hardcoded
  extra_sp = names(pacpal[is.na(pacpal)])
  pacpal = pacpal[c(names(palcols), extra_sp[order(extra_sp)])]

  pacpal[is.na(pacpal)] = colorspace::sequential_hcl(sum(is.na(pacpal)), palette = "Light Grays")

  return(pacpal)
}
