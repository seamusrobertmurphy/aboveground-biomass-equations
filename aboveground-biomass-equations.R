
# Inventory Monitoring and Assessment #
#########  Tree Allometry #############
####### S.Murphy -08-03-2020 ##########

## Load data ##
## Compute HEF per seperately & calculate relascope plots using a:B ratio of 1:35 ##
ForestInventoryFixedSize_Rversion$hef <- 10000 / (pi*ForestInventoryFixedSize_Rversion$plot_radius_m^2)
ForestInventoryConcentric_Rversion$hef <- 10000 / (pi*ForestInventoryConcentric_Rversion$plot_radius_m^2)
ForestInventoryKtree_Rversion$hef <- 10000 / (pi*ForestInventoryKtree_Rversion$plot_radius_m^2)
ForestInventoryRelascope_Rversion$plot_radius_m <- (35*ForestInventoryRelascope_Rversion$dbh_cm/100)
ForestInventoryRelascope_Rversion$hef <- 10000 / (pi*ForestInventoryRelascope_Rversion$plot_radius_m^2)

## Merge datasets ##
total_trees <- rbind(ForestInventoryFixedSize_Rversion, ForestInventoryConcentric_Rversion, ForestInventoryKtree_Rversion, ForestInventoryRelascope_Rversion)

## Generate basal area per tree using dbh ##
total_trees$tree_ba_m <- ((total_trees$dbh_cm)/200)^2*3.147

## split data subsets ##
total_trees$species <- as.factor(total_trees$species)
total_trees$species_code <- recode(total_trees$species, 'SS'=1, 'PBI'=2, 'SBI'=3, 'ROW'=4)
SS_species <- total_trees[ which(total_trees$species_code==1),]
PBI_species <- total_trees[ which(total_trees$species_code==2),]
SBI_species <- total_trees[ which(total_trees$species_code==3),]
ROW_species <- total_trees[ which(total_trees$species_code==4),]

## Generate "agb" response variables using selected biomass equations ##
## Equations used: zianis05.32, zianis05.52, zianis05.53, zianis05.54
## zianis05.55, zianis05.54, bunce68.01, zianis05.40, chojnacky.C14.01.172
## chojnacky.C14.02.172, black.B04.01, chojnacky.C14.16 - see bilbiography ##
SBI_species$tree_agb_kg_Z05.32 <- (((SBI_species$dbh_cm)^2.29)*0.2511)
SBI_species$tree_agb_kg_Z05.52 <- (exp(-2.4166 +2.4227*log(SBI_species$dbh_cm)))
SBI_species$tree_agb_kg_Z05.53 <- (exp(-2.7584 +2.6134*log(SBI_species$dbh_cm)))
SBI_species$tree_agb_kg_Z05.54 <- (exp(-2.1625 +2.3078*log(SBI_species$dbh_cm)))
SBI_species$tree_agb_kg_Z05.55 <- (exp(-2.6423 +2.4678*log(SBI_species$dbh_cm)))
PBI_species$tree_agb_kg_B68.01 <- (exp(-2.162 + 2.3078*log(PBI_species$dbh_cm)))
PBI_species$tree_agb_kg_Z05.40 <- (((PBI_species$dbh_cm)^2.50038)*0.00029)
ROW_species$tree_agb_kg_C14.01 <- (exp(-2.9255 + 2.4109*log(ROW_species$dbh_cm)))
ROW_species$tree_agb_kg_C14.02 <- (exp(-2.2118 + 2.4133*log(ROW_species$dbh_cm)))
SS_species$tree_agb_kg_B04.01 <- (((SS_species$dbh_cm)^2.71)*0.028)
SS_species$tree_agb_kg_C14.16 <- (exp(-3.030 + 2.5567*log(SS_species$dbh_cm)))

## Analysis; Run descriptives and check distrubtion of explanatory/response variables ##
describe(SBI_species$dbh_cm)
describe(PBI_species$dbh_cm)
describe(SS_species$dbh_cm)
describe(ROW_species$dbh_cm)

describe(SBI_species$tree_ba_m)
describe(PBI_species$tree_ba_m)
describe(SS_species$tree_ba_m)
describe(ROW_species$tree_ba_m)

describe(SBI_species$tree_agb_kg_Z05.32) 
describe(SBI_species$tree_agb_kg_Z05.52)
describe(SBI_species$tree_agb_kg_Z05.53)
describe(SBI_species$tree_agb_kg_Z05.54)
describe(SBI_species$tree_agb_kg_Z05.55)
describe(PBI_species$tree_agb_kg_B68.01)
describe(PBI_species$tree_agb_kg_Z05.40)
describe(ROW_species$tree_agb_kg_C14.01)
describe(ROW_species$tree_agb_kg_C14.02)
describe(SS_species$tree_agb_kg_B04.01)
describe(SS_species$tree_agb_kg_C14.16)

## Analysis; Check models for distribution, homoskedasticity##
## Compare linear/non-linear models ##
## Eq z05.32 ##

predict1.z05.32 <- lm(SBI_species$tree_agb_kg_Z05.32 ~ SBI_species$dbh_cm)
plot(SBI_species$tree_agb_kg_Z05.32, resid(predict1.z05.32))
plot(fitted(predict1.z05.32), resid(predict1.z05.32))
qqnorm(rstandard(predict1.z05.32))

plot(SBI_species$dbh_cm, SBI_species$tree_agb_kg_Z05.32)
abline(predict1.z05.32, lwd=3, col="red")
predict2.z05.32 <- lm(SBI_species$tree_agb_kg_Z05.32 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict2.z05.32)), lwd=3, col="blue")
predict3.z05.32 <- lm(SBI_species$tree_agb_kg_Z05.32 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2) + I(SBI_species$dbh_cm^3))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict3.z05.32)), col="green", lwd=3, lty=3)

## Eq z05.52 ##
predict1.z05.52 <- lm(SBI_species$tree_agb_kg_Z05.52 ~ SBI_species$dbh_cm)
plot(SBI_species$tree_agb_kg_Z05.52, resid(predict1.z05.52))
plot(fitted(predict1.z05.52), resid(predict1.z05.52))
qqnorm(rstandard(predict1.z05.52))

plot(SBI_species$dbh_cm, SBI_species$tree_agb_kg_Z05.52)
abline(predict1.z05.52, lwd=3, col="red")
predict2.z05.52 <- lm(SBI_species$tree_agb_kg_Z05.52 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict2.z05.52)), lwd=3, col="blue")
predict3.z05.52 <- lm(SBI_species$tree_agb_kg_Z05.52 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2) + I(SBI_species$dbh_cm^3))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict3.z05.52)), col="green", lwd=3, lty=3)

## Eq z05.53 ##
predict1.z05.53 <- lm(SBI_species$tree_agb_kg_Z05.53 ~ SBI_species$dbh_cm)
plot(SBI_species$tree_agb_kg_Z05.53, resid(predict1.z05.53))
plot(fitted(predict1.z05.53), resid(predict1.z05.53))
qqnorm(rstandard(predict1.z05.53))

plot(SBI_species$dbh_cm, SBI_species$tree_agb_kg_Z05.53)
abline(predict1.z05.53, col="red", lwd=3)
predict2.z05.53 <- lm(SBI_species$tree_agb_kg_Z05.53 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict2.z05.53)), lwd=3, col="blue")
predict3.z05.53 <- lm(SBI_species$tree_agb_kg_Z05.53 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2) + I(SBI_species$dbh_cm^3))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict3.z05.53)), col="green", lwd=3, lty=3)

## Eq z05.54 ##
predict1.z05.54 <- lm(SBI_species$tree_agb_kg_Z05.54 ~ SBI_species$dbh_cm)
plot(SBI_species$tree_agb_kg_Z05.54, resid(predict1.z05.54))
plot(fitted(predict1.z05.54), resid(predict1.z05.54))
qqnorm(rstandard(predict1.z05.54))

plot(SBI_species$dbh_cm, SBI_species$tree_agb_kg_Z05.54)
abline(predict1.z05.54, col="red", lwd=3)
predict2.z05.54 <- lm(SBI_species$tree_agb_kg_Z05.54 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict2.z05.54)), lwd=3, col="blue")
predict3.z05.54 <- lm(SBI_species$tree_agb_kg_Z05.54 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2) + I(SBI_species$dbh_cm^3))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict3.z05.54)), col="green", lwd=3, lty=3)

## Eq z05.55 ##
predict1.z05.55 <- lm(SBI_species$tree_agb_kg_Z05.55 ~ SBI_species$dbh_cm)
plot(SBI_species$tree_agb_kg_Z05.55, resid(predict1.z05.55))
plot(fitted(predict1.z05.55), resid(predict1.z05.55))
qqnorm(rstandard(predict1.z05.55))

plot(SBI_species$dbh_cm, SBI_species$tree_agb_kg_Z05.55)
abline(predict1.z05.55, col="red", lwd=3)
predict2.z05.55 <- lm(SBI_species$tree_agb_kg_Z05.55 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict2.z05.55)), lwd=3, col="blue")
predict3.z05.55 <- lm(SBI_species$tree_agb_kg_Z05.55 ~ SBI_species$dbh_cm + I(SBI_species$dbh_cm^2) + I(SBI_species$dbh_cm^3))
lines(smooth.spline(SBI_species$dbh_cm, predict(predict3.z05.55)), col="green", lwd=3, lty=3)

## Eq B68.01 ##
predict1.B68.01 <- lm(PBI_species$tree_agb_kg_B68.01 ~ PBI_species$dbh_cm)
plot(PBI_species$tree_agb_kg_B68.01, resid(predict1.B68.01))
plot(fitted(predict1.B68.01), resid(predict1.B68.01))
qqnorm(rstandard(predict.B68.01))

plot(PBI_species$dbh_cm, PBI_species$tree_agb_kg_B68.01)
abline(predict1.B68.01, col="red", lwd=3)
predict2.B68.01 <- lm(PBI_species$tree_agb_kg_B68.01 ~ PBI_species$dbh_cm + I(PBI_species$dbh_cm^2))
lines(smooth.spline(PBI_species$dbh_cm, predict(predict2.B68.01)), lwd=3, col="blue")
predict3.B68.01 <- lm(PBI_species$tree_agb_kg_B68.01 ~ PBI_species$dbh_cm + I(PBI_species$dbh_cm^2) + I(PBI_species$dbh_cm^3))
lines(smooth.spline(PBI_species$dbh_cm, predict(predict3.B68.01)), col="green", lwd=3, lty=3)

## Eq Z05.40 ##
predict1.z05.40 <- lm(PBI_species$tree_agb_kg_Z05.40 ~ PBI_species$dbh_cm)
plot(PBI_species$tree_agb_kg_Z05.40, resid(predict1.z05.40))
plot(fitted(predict1.z05.40), resid(predict1.z05.40))
qqnorm(rstandard(predict1.z05.40))

plot(PBI_species$dbh_cm, PBI_species$tree_agb_kg_Z05.40)
abline(predict1.z05.40, col="red", lwd=3)
predict2.z05.40 <- lm(PBI_species$tree_agb_kg_Z05.40 ~ PBI_species$dbh_cm + I(PBI_species$dbh_cm^2))
lines(smooth.spline(PBI_species$dbh_cm, predict(predict2.z05.40)), lwd=3, col="blue")
predict3.z05.40 <- lm(PBI_species$tree_agb_kg_Z05.40 ~ PBI_species$dbh_cm + I(PBI_species$dbh_cm^2) + I(PBI_species$dbh_cm^3))
lines(smooth.spline(PBI_species$dbh_cm, predict(predict3.z05.40)), col="green", lwd=3, lty=3)

## Eq B04.01 ##
predict1.B04.01 <- lm(SS_species$tree_agb_kg_B04.01 ~ SS_species$dbh_cm)
plot(SS_species$tree_agb_kg_B04.01, resid(predict1.B04.01))
plot(fitted(predict1.B04.01), resid(predict1.B04.01))
qqnorm(rstandard(predict1.B04.01))

plot(SS_species$dbh_cm, SS_species$tree_agb_kg_B04.01)
abline(predict1.B04.01, col="red", lwd=3)
predict2.B04.01 <- lm(SS_species$tree_agb_kg_B04.01 ~ SS_species$dbh_cm + I(SS_species$dbh_cm^2))
lines(smooth.spline(SS_species$dbh_cm, predict(predict2.B04.01)), lwd=3, col="blue")
predict3.B04.01 <- lm(SS_species$tree_agb_kg_B04.01 ~ SS_species$dbh_cm + I(SS_species$dbh_cm^2) + I(SS_species$dbh_cm^3))
lines(smooth.spline(SS_species$dbh_cm, predict(predict3.B04.01)), col="green", lwd=3, lty=3)

## Eq C14.16 ##
predict1_C14.16 <- lm(SS_species$tree_agb_kg_C14.16 ~ SS_species$dbh_cm)
plot(SS_species$tree_agb_kg_C14.16, resid(predict1_C14.16))
plot(fitted(predict1_C14.16), resid(predict1_C14.16))
qqnorm(rstandard(predict1_C14.16))

plot(SS_species$dbh_cm, SS_species$tree_agb_kg_C14.16)
abline(predict1_C14.16, col="red", lwd=3)
predict2_C14.16 <- lm(SS_species$tree_agb_kg_C14.16 ~ SS_species$dbh_cm + I(SS_species$dbh_cm^2))
lines(smooth.spline(SS_species$dbh_cm, predict(predict2_C14.16)), lwd=3, col="blue")
predict3_C14.16 <- lm(SS_species$tree_agb_kg_C14.16 ~ SS_species$dbh_cm + I(SS_species$dbh_cm^2) + I(SS_species$dbh_cm^3))
lines(smooth.spline(SS_species$dbh_cm, predict(predict3_C14.16)), col="green", lwd=3, lty=3)

## Compare model resuls ## 
summary(predict2.z05.32)
summary(predict2.z05.52)
summary(predict2.z05.53)
summary(predict2.z05.54)
summary(predict2.z05.55)
summary(predict2.B68.01)
summary(predict2.z05.40)
summary(predict2.B04.01)
summary(predict2_C14.16)

## split fixed_size_data to estiamte plot level BA & AGB##
sapply(ForestInventoryFixedSize_Rversion, class)
ForestInventoryFixedSize_Rversion$species <- as.factor(ForestInventoryFixedSize_Rversion$species)
ForestInventoryFixedSize_Rversion$species_code <- recode(ForestInventoryFixedSize_Rversion$species, 'SS'=1, 'PBI'=2, 'SBI'=3, 'ROW'=4)
SS_species <- ForestInventoryFixedSize_Rversion[ which(ForestInventoryFixedSize_Rversion$species_code==1),]
PBI_species <- ForestInventoryFixedSize_Rversion[ which(ForestInventoryFixedSize_Rversion$species_code==2),]
SBI_species <- ForestInventoryFixedSize_Rversion[ which(ForestInventoryFixedSize_Rversion$species_code==3),]
ROW_species <- ForestInventoryFixedSize_Rversion[ which(ForestInventoryFixedSize_Rversion$species_code==4),]

## generate tree_agb_kg using selected biomass equations ##
SBI_species$tree_agb_kg <- (exp(-2.7584 +2.6134*log(SBI_species$dbh_cm)))
PBI_species$tree_agb_kg <- (((PBI_species$dbh_cm)^2.50038)*0.00029)
ROW_species$tree_agb_kg <- (exp(-2.9255 + 2.4109*log(ROW_species$dbh_cm)))
SS_species$tree_agb_kg <- (exp(-3.030 + 2.5567*log(SS_species$dbh_cm)))

# merge fixed_size_datasets 
fixedplot_merged <- rbind(PBI_species, SBI_species, ROW_species, SS_species)

fixedplot_merged$tree_ba_m <- ((fixedplot_merged$dbh_cm)/200)^2*3.147
fixedplot_merged$BA_m <- (fixedplot_merged$tree_ba_m*fixedplot_merged$hef)
fixedplot_merged$AGB_Mg_ha <- (fixedplot_merged$tree_agb_kg*fixedplot_merged$hef)

fixedplot_merged %>%
  group_by(plot) %>%
  summarise(sumplot_AGB = sum(fixedplot_merged$BA_m, na.rm = TRUE))

