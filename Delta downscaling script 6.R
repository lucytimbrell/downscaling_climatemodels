## Timbrell et al. (2022) Comparison between climate model and proxy reconstructions for Late Pleistocene climate

################# PREPARATION  ################

##  Install pastlclim and terra

devtools::install_github("EvolEcolGroup/pastclim",ref="ds_files") #  CHECK WITH ANDREA WHICH ONE TO INCLUDE IN FINAL SCRIPT
install.packages('terra', repos='https://rspatial.r-universe.dev')
remotes::install_github("ericpante/marmap")

##  Load packages 
library(terra) 
library(pastclim)

##  Set working directory 
setwd("~/Documents/LT_documents/publications/biology letters")

##  Set spatial extent 
sea_ext<- terra::ext(-180, 180, 0, 90) # Northern hempisphere

################# LOAD DATA  ################
## Load in reconstructions 

download_dataset(data="Beyer2020") # download Beyer 2020

## Proxy data

proxies <- read.csv("Final_proxy_data2.csv")
proxies$Age..ka.BP...mean. <-  proxies$Age..ka.BP...mean.*1000
proxies$Age..ka.BP...mean. <-  proxies$Age..ka.BP...mean.*(-1) # Negative 
proxies$Age..ka.BP...mean.rounded<- round(proxies$Age..ka.BP...mean., digits = -3) # Round to nearest 1000 to match chronological resolution

# Time steps needed for downscaling

time_steps <- rev(get_time_steps(dataset = "Beyer2020")) # Check if timeslices are available

################# DELTA DOWNSCALING - TEMPERATURE ################

tavg_vars <- c(paste0("temperature_0",1:9),paste0("temperature_",10:12))

tavg_series <- region_series(time_bp = time_steps, 
                             bio_variables = tavg_vars, 
                             dataset = 'Beyer2020', 
                             ext = sea_ext) 

## Get temperature to compute high res masks
tavg_model_lres_rast <- tavg_series$temperature_01
tavg_model_lres_rast

##  Get high res modern data
pastclim:::download_worldclim("tavg",10)

tavg_obs_hres_all<- pastclim:::load_worldclim("tavg",10)
tavg_obs_hres_all <- terra::crop(tavg_obs_hres_all, sea_ext)

## Get first month
tavg_obs_hres_rast <- tavg_obs_hres_all[[1]]

ext(tavg_obs_hres_rast)==ext(tavg_model_lres_rast)## check extents


##  Get relief raster

# Import bathymetric map for reconstruction of sea level change

relief_rast <- pastclim:::download_relief(tavg_obs_hres_rast)

ext(relief_rast)==ext(tavg_model_lres_rast ) # Check same extent
res(tavg_obs_hres_rast) == res(relief_rast) # Check same resolution
ncol(relief_rast) == ncol(tavg_obs_hres_rast)
nrow(relief_rast) == nrow(tavg_obs_hres_rast)


# Make land masks for each time step
high_res_mask <- pastclim:::make_land_mask(relief_rast = relief_rast,
                                 time_bp = time_bp(tavg_model_lres_rast))
plot(high_res_mask[[1:12]], col = "blue") # Plot to check 

ext(high_res_mask)==ext(tavg_obs_hres_all) # Check same extent

## Delta-downscaling for all months 
tavg_downscaled_list<-list()
for (i in 1:12){
  delta_rast<-pastclim:::delta_compute(x=tavg_series[[i]], ref_time = 0, 
                                       obs = tavg_obs_hres_all[[i]])
  tavg_downscaled_list[[i]] <- pastclim:::delta_downscale (x = tavg_series[[i]], 
                                                           delta_rast = delta_rast,
                                                           x_landmask_high = high_res_mask)
  }
tavg_downscaled <- terra::sds(tavg_downscaled_list)

tavg_downscaled # inspect

################# DELTA DOWNSCALING - PRECIPITATION  ################
prec_vars <- c(paste0("precipitation_0",1:9),paste0("precipitation_",10:12))
prec_series <- region_series(bio_variables = prec_vars, # Extract regional time series 
                             time_bp =  time_steps,
                             dataset = "Beyer2020",
                             ext = sea_ext)

pastclim:::download_worldclim("prec",10) #  High res reconstructions
prec_obs_hres_all <- pastclim:::load_worldclim("prec",10)
prec_obs_hres_all <- terra::crop(prec_obs_hres_all, sea_ext) #  Crop

## Delta downscaling
prec_downscaled_list<-list() 

for (i in 1:12){
  delta_rast<-pastclim:::delta_compute(x=prec_series[[i]], ref_time = 0,
                                       obs = prec_obs_hres_all[[i]])
  prec_downscaled_list[[i]] <- pastclim:::delta_downscale (x = prec_series[[i]],
                                                            delta_rast = delta_rast,
                                                           x_landmask_high = high_res_mask)
print(i)
  }
prec_downscaled <- terra::sds(prec_downscaled_list)

################# COMPUTE BIOVARIABLES  ################
bioclim_downscaled<-bioclim_vars(tavg = tavg_downscaled, prec = prec_downscaled)

plot(bioclim_downscaled[[1]], main = time(bioclim_downscaled[[1]]))

################# SAVE DATASETS ################

setwd() # SET TO HARD DRIVE 

terra::writeCDF(bioclim_downscaled[[1]], "bio01_downscaled.nc", varname = "bio01", overwrite = TRUE) # mean annual temp
terra::writeCDF(bioclim_downscaled[[2]], "bio04_downscaled.nc", varname = "bio04", overwrite = TRUE) # temp seasonality 
terra::writeCDF(bioclim_downscaled[[3]], "bio05_downscaled.nc", varname = "bio05", overwrite = TRUE) # min annual temp
terra::writeCDF(bioclim_downscaled[[4]], "bio06_downscaled.nc", varname = "bio06", overwrite = TRUE) # max annual temp
terra::writeCDF(bioclim_downscaled[[5]], "bio07_downscaled.nc", varname = "bio07", overwrite = TRUE) # temp annual range
terra::writeCDF(bioclim_downscaled[[6]], "bio08_downscaled.nc", varname = "bio08", overwrite = TRUE) # mean temp of wettest quarter
terra::writeCDF(bioclim_downscaled[[7]], "bio09_downscaled.nc", varname = "bio09", overwrite = TRUE) # mean temp of driest quarter
terra::writeCDF(bioclim_downscaled[[8]], "bio10_downscaled.nc", varname = "bio10", overwrite = TRUE) # mean temp of warmest quarter
terra::writeCDF(bioclim_downscaled[[9]], "bio11_downscaled.nc", varname = "bio11", overwrite = TRUE) # mean temp of coldest quarter
terra::writeCDF(bioclim_downscaled[[10]], "bio12_downscaled.nc", varname = "bio12", overwrite = TRUE) # total annual precip
terra::writeCDF(bioclim_downscaled[[11]], "bio13_downscaled.nc", varname = "bio13", overwrite = TRUE) # precip of wettest month
terra::writeCDF(bioclim_downscaled[[12]], "bio14_downscaled.nc", varname = "bio14", overwrite = TRUE) # precip of driest month
terra::writeCDF(bioclim_downscaled[[13]], "bio15_downscaled.nc", varname = "bio15", overwrite = TRUE) # precip seasonality
terra::writeCDF(bioclim_downscaled[[14]], "bio16_downscaled.nc", varname = "bio16", overwrite = TRUE) # precip of wettest quarter
terra::writeCDF(bioclim_downscaled[[15]], "bio17_downscaled.nc", varname = "bio17", overwrite = TRUE) # precip of driest quarter
terra::writeCDF(bioclim_downscaled[[16]], "bio18_downscaled.nc", varname = "bio18", overwrite = TRUE) # precip of warmest quarter
terra::writeCDF(bioclim_downscaled[[17]], "bio19_downscaled.nc", varname = "bio19", overwrite = TRUE) #  precip of coldest quarter

################# RESET TIME UNITS ################

nc_name <- "bio01_downscaled.nc"
nc_in <- ncdf4::nc_open(nc_name, write=TRUE)
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "units",
                 attval = "years since 1950-01-01 00:00:00.0")
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "long_name",
                 attval = "years BP")
ncdf4::ncatt_put(nc_in, varid="time", attname="axis", attval = "T")
ncdf4::nc_close(nc_in)

nc_name <- "bio12_downscaled.nc"
nc_in <- ncdf4::nc_open(nc_name, write=TRUE)
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "units",
                 attval = "years since 1950-01-01 00:00:00.0")
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "long_name",
                 attval = "years BP")
ncdf4::ncatt_put(nc_in, varid="time", attname="axis", attval = "T")
ncdf4::nc_close(nc_in)

nc_name <- "bio10_downscaled.nc"
nc_in <- ncdf4::nc_open(nc_name, write=TRUE)
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "units",
                 attval = "years since 1950-01-01 00:00:00.0")
ncdf4::ncatt_put(nc_in,varid="time", 
                 attname = "long_name",
                 attval = "years BP")
ncdf4::ncatt_put(nc_in, varid="time", attname="axis", attval = "T")
ncdf4::nc_close(nc_in)

################# EXTRACT TIME SERIES AND PLOT ################
library(ggplot2)
library(ggpubr)

## Proxy 1 - Luanhaizi LH2
coords <- data.frame(longitude=101.35,latitude=c(37.593333))
proxy <- subset(proxies, Site == "Luanhaizi LH2")
time_steps <- unique(sort(proxy$Age..ka.BP...mean.rounded, decreasing = TRUE))
time_steps <- time_steps[which(time_steps %in% get_time_steps(dataset = "Beyer2020"))] # Check if timeslices are available

  
raw_time_series_temp <-location_series(x=coords,time_bp=time_steps, # raw data
                                   bio_variables=c("bio01"),
                                   dataset="Beyer2020")

downscaled_time_series_temp <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                           bio_variables=c("bio01"),
                                           dataset="custom", path_to_nc = "bio01_downscaled.nc")

raw_time_series_prec <-location_series(x=coords,time_bp=time_steps, # raw data
                                        bio_variables=c("bio12"),
                                        dataset="Beyer2020")

downscaled_time_series_prec <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                                bio_variables=c("bio12"),
                                                dataset="custom", path_to_nc = "bio12_downscaled.nc")


a <- ggplot() + 
  geom_line(data = raw_time_series_temp, aes(x= time, y = bio01, colour = "blue")) + 
  geom_line(data = downscaled_time_series_temp, aes(x = time, y = bio01, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Luanhaizi LH2"), aes(x=Age..ka.BP...mean., y = Annual.WA.PLS, colour = "red"))+
  geom_line(data =subset(proxies, Site == "Luanhaizi LH2"), aes(x=Age..ka.BP...mean., y = Annual.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Temperature (deg. c)") +
  ggtitle("Luanhaizi LH2")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")
  

b <-  ggplot() + 
  geom_line(data = raw_time_series_prec, aes(x= time, y = bio12, colour = "blue")) + 
  geom_line(data = downscaled_time_series_prec, aes(x = time, y = bio12, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Luanhaizi LH2"), aes(x=Age..ka.BP...mean., y = Total.precip.WA.PLS, colour = "red"))+
  geom_line(data = subset(proxies, Site == "Luanhaizi LH2"), aes(x=Age..ka.BP...mean., y = Total.precip.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Precipitation (mm)") +
  #ggtitle("Luanhaizi LH2")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")

#ggarrange(a, b, nrow = 1, top = "Luanhaizi LH2")

## Calculate bias and RMSE between time series 

Time <- raw_time_series_temp$time
 
new_MAT_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.MAT, xout = Time, method = "linear") # interpolate
new_WA.PLS_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.WA.PLS, xout = Time, method = "linear")
new_MAT_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.MAT, xout = Time, method = "linear")
new_WA.PLS_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.WA.PLS, xout = Time, method = "linear")

df <- data.frame(Time = Time, 
                 RawT = raw_time_series_temp$bio01, 
                 DownscaledT = downscaled_time_series_temp$bio01,
                 MAT_temp = new_MAT_temp$y,
                 WA.PLS_temp = new_WA.PLS_temp$y,
                 RawP = raw_time_series_prec$bio12, 
                 DownscaledP = downscaled_time_series_prec$bio12,
                 MAT_prec = new_MAT_prec$y,
                 WA.PLS_prec = new_WA.PLS_prec$y)

proxy1_rawT_MAT_bias <- mean(df$RawT - df$MAT_temp)
proxy1_rawT_MAT_RMSE <- sqrt(mean((df$RawT - df$MAT_temp)^2))
proxy1_rawT_WAPLS_bias <- mean(df$RawT - df$WA.PLS_temp)
proxy_1rawT_WAPLS_RMSE <- sqrt(mean((df$RawT - df$WA.PLS_temp)^2))

proxy1_dsT_MAT_bias <- mean(df$DownscaledT - df$MAT_temp)
proxy1_dsT_MAT_RMSE <- sqrt(mean((df$DownscaledT - df$MAT_temp)^2))
proxy1_dsT_WAPLS_bias <- mean(df$DownscaledT - df$WA.PLS_temp)
proxy_1dsT_WAPLS_RMSE <- sqrt(mean((df$DownscaledT - df$WA.PLS_temp)^2))

proxy1_rawP_MAT_bias <- mean(df$RawP - df$MAT_prec)
proxy1_rawP_MAT_RMSE <- sqrt(mean((df$RawP - df$MAT_prec)^2))
proxy1_rawP_WAPLS_bias <- mean(df$RawP - df$WA.PLS_prec)
proxy_1rawP_WAPLS_RMSE <- sqrt(mean((df$RawP - df$WA.PLS_prec)^2))

proxy1_dsP_MAT_bias <- mean(df$DownscaledP - df$MAT_prec)
proxy1_dsP_MAT_RMSE <- sqrt(mean((df$DownscaledP - df$MAT_prec)^2))
proxy1_dsP_WAPLS_bias <- mean(df$DownscaledP - df$WA.PLS_prec)
proxy_1dsP_WAPLS_RMSE <- sqrt(mean((df$DownscaledP - df$WA.PLS_prec)^2))


## Proxy 2 - Pittsburg Basin
coords <- data.frame(longitude=-89.18984,latitude=c(38.90363))
proxy <- subset(proxies, Site == "Pittsburg Basin")
time_steps <- unique(sort(proxy$Age..ka.BP...mean.rounded, decreasing = TRUE))
time_steps <- time_steps[which(time_steps %in% get_time_steps(dataset = "Beyer2020"))] # Check if timeslices are available

  
raw_time_series_temp <-location_series(x=coords,time_bp=time_steps, # raw data
                                         bio_variables=c("bio01"),
                                         dataset="Beyer2020")

downscaled_time_series_temp <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio01"),
                                               dataset="custom", path_to_nc = "bio01_downscaled.nc")

raw_time_series_prec <-location_series(x=coords,time_bp=time_steps, # raw data
                                       bio_variables=c("bio12"),
                                       dataset="Beyer2020")

downscaled_time_series_prec <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio12"),
                                               dataset="custom", path_to_nc = "bio12_downscaled.nc")

c <- ggplot() + 
  geom_line(data = raw_time_series_temp, aes(x= time, y = bio01, colour = "blue")) + 
  geom_line(data = downscaled_time_series_temp, aes(x = time, y = bio01, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Pittsburg Basin"), aes(x=Age..ka.BP...mean., y = Annual.WA.PLS, colour = "red"))+
  geom_line(data =subset(proxies, Site == "Pittsburg Basin"), aes(x=Age..ka.BP...mean., y = Annual.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Temperature (deg. c)") +
  ggtitle("Pittsburg Basin")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")

d <-  ggplot() + 
  geom_line(data = raw_time_series_prec, aes(x= time, y = bio12, colour = "blue")) + 
  geom_line(data = downscaled_time_series_prec, aes(x = time, y = bio12, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Pittsburg Basin"), aes(x=Age..ka.BP...mean., y = Total.precip.WA.PLS, colour = "red"))+
  geom_line(data = subset(proxies, Site == "Pittsburg Basin"), aes(x=Age..ka.BP...mean., y = Total.precip.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Precipitation (mm)") +
  #ggtitle("Pittsburg Basin")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")

#ggarrange(c, d, nrow = 1, top = "Pittsburg Basin")

## Calculate bias and RMSE between time series 

Time <- raw_time_series_temp$time

new_MAT_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.MAT, xout = Time, method = "linear") # interpolate
new_WA.PLS_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.WA.PLS, xout = Time, method = "linear")
new_MAT_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.MAT, xout = Time, method = "linear")
new_WA.PLS_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.WA.PLS, xout = Time, method = "linear")

df <- data.frame(Time = Time, 
                 RawT = raw_time_series_temp$bio01, 
                 DownscaledT = downscaled_time_series_temp$bio01,
                 MAT_temp = new_MAT_temp$y,
                 WA.PLS_temp = new_WA.PLS_temp$y,
                 RawP = raw_time_series_prec$bio12, 
                 DownscaledP = downscaled_time_series_prec$bio12,
                 MAT_prec = new_MAT_prec$y,
                 WA.PLS_prec = new_WA.PLS_prec$y)

proxy2_rawT_MAT_bias <- mean(df$RawT - df$MAT_temp, na.rm = TRUE)
proxy2_rawT_MAT_RMSE <- sqrt(mean((df$RawT - df$MAT_temp)^2, na.rm = TRUE))
proxy2_rawT_WAPLS_bias <- mean(df$RawT - df$WA.PLS_temp, na.rm = TRUE)
proxy_2rawT_WAPLS_RMSE <- sqrt(mean((df$RawT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy2_dsT_MAT_bias <- mean(df$DownscaledT - df$MAT_temp, na.rm = TRUE)
proxy2_dsT_MAT_RMSE <- sqrt(mean((df$DownscaledT - df$MAT_temp)^2, na.rm = TRUE))
proxy2_dsT_WAPLS_bias <- mean(df$DownscaledT - df$WA.PLS_temp, na.rm = TRUE)
proxy2_dsT_WAPLS_RMSE <- sqrt(mean((df$DownscaledT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy2_rawP_MAT_bias <- mean(df$RawP - df$MAT_prec, na.rm = TRUE)
proxy2_rawP_MAT_RMSE <- sqrt(mean((df$RawP - df$MAT_prec)^2, na.rm = TRUE))
proxy2_rawP_WAPLS_bias <- mean(df$RawP - df$WA.PLS_prec, na.rm = TRUE)
proxy2_rawP_WAPLS_RMSE <- sqrt(mean((df$RawP - df$WA.PLS_prec)^2, na.rm = TRUE))

proxy2_dsP_MAT_bias <- mean(df$DownscaledP - df$MAT_prec, na.rm = TRUE)
proxy2_dsP_MAT_RMSE <- sqrt(mean((df$DownscaledP - df$MAT_prec)^2, na.rm = TRUE))
proxy2_dsP_WAPLS_bias <- mean(df$DownscaledP - df$WA.PLS_prec, na.rm = TRUE)
proxy2_dsP_WAPLS_RMSE <- sqrt(mean((df$DownscaledP - df$WA.PLS_prec)^2, na.rm = TRUE))

## Proxy 3 - Hanging Lake
coords <- data.frame(longitude=-138.383333,latitude=c(68.383333))
proxy <- subset(proxies, Site == "Hanging Lake")
time_steps <- unique(sort(proxy$Age..ka.BP...mean.rounded, decreasing = TRUE))
time_steps <- time_steps[which(time_steps %in% get_time_steps(dataset = "Beyer2020"))] # Check if timeslices are available

  raw_time_series_temp <-location_series(x=coords,time_bp=time_steps, # raw data
                                         bio_variables=c("bio01"),
                                         dataset="Beyer2020")

downscaled_time_series_temp <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio01"),
                                               dataset="custom", path_to_nc = "bio01_downscaled.nc")

raw_time_series_prec <-location_series(x=coords,time_bp=time_steps, # raw data
                                       bio_variables=c("bio12"),
                                       dataset="Beyer2020")

downscaled_time_series_prec <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio12"),
                                               dataset="custom", path_to_nc = "bio12_downscaled.nc")

e <- ggplot() + 
  geom_line(data = raw_time_series_temp, aes(x= time, y = bio01, colour = "blue")) + 
  geom_line(data = downscaled_time_series_temp, aes(x = time, y = bio01, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Hanging Lake"), aes(x=Age..ka.BP...mean., y = Annual.WA.PLS, colour = "red"))+
  geom_line(data =subset(proxies, Site == "Hanging Lake"), aes(x=Age..ka.BP...mean., y = Annual.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Temperature (deg. c)") +
  ggtitle("Hanging Lake")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")



f <-  ggplot() + 
  geom_line(data = raw_time_series_prec, aes(x= time, y = bio12, colour = "blue")) + 
  geom_line(data = downscaled_time_series_prec, aes(x = time, y = bio12, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Hanging Lake"), aes(x=Age..ka.BP...mean., y = Total.precip.WA.PLS, colour = "red"))+
  geom_line(data = subset(proxies, Site == "Hanging Lake"), aes(x=Age..ka.BP...mean., y = Total.precip.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
   ylab("Precipitation (mm)") +
  # ggtitle("Hanging Lake")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")

#grid.arrange(e, f, nrow = 1, top = "Hanging Lake")


## Calculate bias and RMSE between time series 

Time <- raw_time_series_temp$time

new_MAT_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.MAT, xout = Time, method = "linear") # interpolate
new_WA.PLS_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.WA.PLS, xout = Time, method = "linear")
new_MAT_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.MAT, xout = Time, method = "linear")
new_WA.PLS_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.WA.PLS, xout = Time, method = "linear")

df <- data.frame(Time = Time, 
                 RawT = raw_time_series_temp$bio01, 
                 DownscaledT = downscaled_time_series_temp$bio01,
                 MAT_temp = new_MAT_temp$y,
                 WA.PLS_temp = new_WA.PLS_temp$y,
                 RawP = raw_time_series_prec$bio12, 
                 DownscaledP = downscaled_time_series_prec$bio12,
                 MAT_prec = new_MAT_prec$y,
                 WA.PLS_prec = new_WA.PLS_prec$y)

proxy3_rawT_MAT_bias <- mean(df$RawT - df$MAT_temp, na.rm = TRUE)
proxy3_rawT_MAT_RMSE <- sqrt(mean((df$RawT - df$MAT_temp)^2, na.rm = TRUE))
proxy3_rawT_WAPLS_bias <- mean(df$RawT - df$WA.PLS_temp, na.rm = TRUE)
proxy3_rawT_WAPLS_RMSE <- sqrt(mean((df$RawT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy3_dsT_MAT_bias <- mean(df$DownscaledT - df$MAT_temp, na.rm = TRUE)
proxy3_dsT_MAT_RMSE <- sqrt(mean((df$DownscaledT - df$MAT_temp)^2, na.rm = TRUE))
proxy3_dsT_WAPLS_bias <- mean(df$DownscaledT - df$WA.PLS_temp, na.rm = TRUE)
proxy3_dsT_WAPLS_RMSE <- sqrt(mean((df$DownscaledT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy3_rawP_MAT_bias <- mean(df$RawP - df$MAT_prec, na.rm = TRUE)
proxy3_rawP_MAT_RMSE <- sqrt(mean((df$RawP - df$MAT_prec)^2, na.rm = TRUE))
proxy3_rawP_WAPLS_bias <- mean(df$RawP - df$WA.PLS_prec, na.rm = TRUE)
proxy3_rawP_WAPLS_RMSE <- sqrt(mean((df$RawP - df$WA.PLS_prec)^2, na.rm = TRUE))

proxy3_dsP_MAT_bias <- mean(df$DownscaledP - df$MAT_prec, na.rm = TRUE)
proxy3_dsP_MAT_RMSE <- sqrt(mean((df$DownscaledP - df$MAT_prec)^2, na.rm = TRUE))
proxy3_dsP_WAPLS_bias <- mean(df$DownscaledP - df$WA.PLS_prec, na.rm = TRUE)
proxy3_dsP_WAPLS_RMSE <- sqrt(mean((df$DownscaledP - df$WA.PLS_prec)^2, na.rm = TRUE))

## Proxy 4 - Elikchan 4 Lake
coords <- data.frame(longitude=151.883333,latitude=c(60.75))
proxy <- subset(proxies, Site == "Elikchan 4 Lake")
time_steps <- unique(sort(proxy$Age..ka.BP...mean.rounded, decreasing = TRUE))
time_steps <- time_steps[which(time_steps %in% get_time_steps(dataset = "Beyer2020"))] # Check if timeslices are available

  raw_time_series_temp <-location_series(x=coords,time_bp=time_steps, # raw data
                                         bio_variables=c("bio01"),
                                         dataset="Beyer2020")

downscaled_time_series_temp <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio01"),
                                               dataset="custom", path_to_nc = "bio01_downscaled.nc")

raw_time_series_prec <-location_series(x=coords,time_bp=time_steps, # raw data
                                       bio_variables=c("bio12"),
                                       dataset="Beyer2020")

downscaled_time_series_prec <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio12"),
                                               dataset="custom", path_to_nc = "bio12_downscaled.nc")

g <- ggplot() + 
  geom_line(data = raw_time_series_temp, aes(x= time, y = bio01, colour = "blue")) + 
  geom_line(data = downscaled_time_series_temp, aes(x = time, y = bio01, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Elikchan 4 Lake"), aes(x=Age..ka.BP...mean., y = Annual.WA.PLS, colour = "red"))+
  geom_line(data =subset(proxies, Site == "Elikchan 4 Lake"), aes(x=Age..ka.BP...mean., y = Annual.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  xlim(min(time_steps), max(time_steps))+
  ylab("Temperature (deg. c)") +
  ggtitle("Elikchan 4 Lake")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")



h <-  ggplot() + 
  geom_line(data = raw_time_series_prec, aes(x= time, y = bio12, colour = "blue")) + 
  geom_line(data = downscaled_time_series_prec, aes(x = time, y = bio12, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Elikchan 4 Lake"), aes(x=Age..ka.BP...mean., y = Total.precip.WA.PLS, colour = "red"))+
  geom_line(data = subset(proxies, Site == "Elikchan 4 Lake"), aes(x=Age..ka.BP...mean., y = Total.precip.MAT, colour = "darkorange"))+
  xlab("year calKa BP") +
  xlim(min(time_steps), max(time_steps))+
  ylab("Precipitation (mm)") +
  # ggtitle("Tianshuihai TS95")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")


ggarrange(g, h, nrow = 1, top = "Elikchan 4 Lake")

## Calculate bias and RMSE between time series 

Time <- raw_time_series_temp$time

new_MAT_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.MAT, xout = Time, method = "linear") # interpolate
new_WA.PLS_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.WA.PLS, xout = Time, method = "linear")
new_MAT_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.MAT, xout = Time, method = "linear")
new_WA.PLS_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.WA.PLS, xout = Time, method = "linear")

df <- data.frame(Time = Time, 
                 RawT = raw_time_series_temp$bio01, 
                 DownscaledT = downscaled_time_series_temp$bio01,
                 MAT_temp = new_MAT_temp$y,
                 WA.PLS_temp = new_WA.PLS_temp$y,
                 RawP = raw_time_series_prec$bio12, 
                 DownscaledP = downscaled_time_series_prec$bio12,
                 MAT_prec = new_MAT_prec$y,
                 WA.PLS_prec = new_WA.PLS_prec$y)

proxy4_rawT_MAT_bias <- mean(df$RawT - df$MAT_temp, na.rm = TRUE)
proxy4_rawT_MAT_RMSE <- sqrt(mean((df$RawT - df$MAT_temp)^2, na.rm = TRUE))
proxy4_rawT_WAPLS_bias <- mean(df$RawT - df$WA.PLS_temp, na.rm = TRUE)
proxy4_rawT_WAPLS_RMSE <- sqrt(mean((df$RawT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy4_dsT_MAT_bias <- mean(df$DownscaledT - df$MAT_temp, na.rm = TRUE)
proxy4_dsT_MAT_RMSE <- sqrt(mean((df$DownscaledT - df$MAT_temp)^2, na.rm = TRUE))
proxy4_dsT_WAPLS_bias <- mean(df$DownscaledT - df$WA.PLS_temp, na.rm = TRUE)
proxy4_dsT_WAPLS_RMSE <- sqrt(mean((df$DownscaledT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy4_rawP_MAT_bias <- mean(df$RawP - df$MAT_prec, na.rm = TRUE)
proxy4_rawP_MAT_RMSE <- sqrt(mean((df$RawP - df$MAT_prec)^2, na.rm = TRUE))
proxy4_rawP_WAPLS_bias <- mean(df$RawP - df$WA.PLS_prec, na.rm = TRUE)
proxy4_rawP_WAPLS_RMSE <- sqrt(mean((df$RawP - df$WA.PLS_prec)^2, na.rm = TRUE))

proxy4_dsP_MAT_bias <- mean(df$DownscaledP - df$MAT_prec, na.rm = TRUE)
proxy4_dsP_MAT_RMSE <- sqrt(mean((df$DownscaledP - df$MAT_prec)^2, na.rm = TRUE))
proxy4_dsP_WAPLS_bias <- mean(df$DownscaledP - df$WA.PLS_prec, na.rm = TRUE)
proxy4_dsP_WAPLS_RMSE <- sqrt(mean((df$DownscaledP - df$WA.PLS_prec)^2, na.rm = TRUE))

## Proxy 5 - Lago Grande di Monticchio
coords <- data.frame(longitude=15.605185,latitude=c(40.931235))
proxy <- subset(proxies, Site == "Lago Grande di Monticchio")
time_steps <- unique(sort(proxy$Age..ka.BP...mean.rounded, decreasing = TRUE))
time_steps <- time_steps[which(time_steps %in% get_time_steps(dataset = "Beyer2020"))] # Check if timeslices are available

  
raw_time_series_temp <-location_series(x=coords,time_bp=time_steps, # raw data
                                         bio_variables=c("bio01"),
                                         dataset="Beyer2020")

downscaled_time_series_temp <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio01"),
                                               dataset="custom", path_to_nc = "bio01_downscaled.nc")

raw_time_series_prec <-location_series(x=coords,time_bp=time_steps, # raw data
                                       bio_variables=c("bio12"),
                                       dataset="Beyer2020")

downscaled_time_series_prec <- location_series(x=coords, time_bp=time_steps, # downscaled data
                                               bio_variables=c("bio12"),
                                               dataset="custom", path_to_nc = "bio12_downscaled.nc")

i <- ggplot() + 
  geom_line(data = raw_time_series_temp, aes(x= time, y = bio01, colour = "blue")) + 
  geom_line(data = downscaled_time_series_temp, aes(x = time, y = bio01, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Lago Grande di Monticchio"), aes(x=Age..ka.BP...mean., y = Annual.WA.PLS, colour = "red"))+
  geom_line(data =subset(proxies, Site == "Lago Grande di Monticchio"), aes(x=Age..ka.BP...mean., y = Annual.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Temperature (deg. c)") +
  ggtitle("Lago Grande di Monticchio")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")


j <-  ggplot() + 
  geom_line(data = raw_time_series_prec, aes(x= time, y = bio12, colour = "blue")) + 
  geom_line(data = downscaled_time_series_prec, aes(x = time, y = bio12, colour = "gray40")) +
  geom_line(data = subset(proxies, Site == "Lago Grande di Monticchio"), aes(x=Age..ka.BP...mean., y = Total.precip.WA.PLS, colour = "red"))+
  geom_line(data = subset(proxies, Site == "Lago Grande di Monticchio"), aes(x=Age..ka.BP...mean., y = Total.precip.MAT, colour = "darkorange"))+
  xlab("year calKa BP") + 
  ylab("Precipitation (mm)") +
  #  ggtitle("Lago Grande di Monticchio")+
  labs(color = "") + 
  scale_color_manual(breaks = c("blue", "gray40", "red", "darkorange"), values = c("blue", "gray40", "red", "darkorange"), labels = c("Raw model data", "Downscaled model data", "Proxy data - WA.PLS method", "Proxy data - MAT method"))+
  theme_bw()+
  theme(legend.position = "none")

ggarrange(a,b,c,d,e,f,g,h,i,j, nrow = 5, ncol = 2, common.legend = TRUE)
ggarrange(e,f,c,d,i,j,a,b,g,h, nrow = 5, ncol = 2, common.legend = TRUE)

## Calculate bias and RMSE between time series 

Time <- raw_time_series_temp$time

new_MAT_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.MAT, xout = Time, method = "linear") # interpolate
new_WA.PLS_temp <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Annual.WA.PLS, xout = Time, method = "linear")
new_MAT_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.MAT, xout = Time, method = "linear")
new_WA.PLS_prec <- approx(x = proxy$Age..ka.BP...mean., y = proxy$Total.precip.WA.PLS, xout = Time, method = "linear")

df <- data.frame(Time = Time, 
                 RawT = raw_time_series_temp$bio01, 
                 DownscaledT = downscaled_time_series_temp$bio01,
                 MAT_temp = new_MAT_temp$y,
                 WA.PLS_temp = new_WA.PLS_temp$y,
                 RawP = raw_time_series_prec$bio12, 
                 DownscaledP = downscaled_time_series_prec$bio12,
                 MAT_prec = new_MAT_prec$y,
                 WA.PLS_prec = new_WA.PLS_prec$y)

proxy5_rawT_MAT_bias <- mean(df$RawT - df$MAT_temp, na.rm = TRUE)
proxy5_rawT_MAT_RMSE <- sqrt(mean((df$RawT - df$MAT_temp)^2, na.rm = TRUE))
proxy5_rawT_WAPLS_bias <- mean(df$RawT - df$WA.PLS_temp, na.rm = TRUE)
proxy5_rawT_WAPLS_RMSE <- sqrt(mean((df$RawT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy5_dsT_MAT_bias <- mean(df$DownscaledT - df$MAT_temp, na.rm = TRUE)
proxy5_dsT_MAT_RMSE <- sqrt(mean((df$DownscaledT - df$MAT_temp)^2, na.rm = TRUE))
proxy5_dsT_WAPLS_bias <- mean(df$DownscaledT - df$WA.PLS_temp, na.rm = TRUE)
proxy5_dsT_WAPLS_RMSE <- sqrt(mean((df$DownscaledT - df$WA.PLS_temp)^2, na.rm = TRUE))

proxy5_rawP_MAT_bias <- mean(df$RawP - df$MAT_prec, na.rm = TRUE)
proxy5_rawP_MAT_RMSE <- sqrt(mean((df$RawP - df$MAT_prec)^2, na.rm = TRUE))
proxy5_rawP_WAPLS_bias <- mean(df$RawP - df$WA.PLS_prec, na.rm = TRUE)
proxy5_rawP_WAPLS_RMSE <- sqrt(mean((df$RawP - df$WA.PLS_prec)^2, na.rm = TRUE))

proxy5_dsP_MAT_bias <- mean(df$DownscaledP - df$MAT_prec, na.rm = TRUE)
proxy5_dsP_MAT_RMSE <- sqrt(mean((df$DownscaledP - df$MAT_prec)^2, na.rm = TRUE))
proxy5_dsP_WAPLS_bias <- mean(df$DownscaledP - df$WA.PLS_prec, na.rm = TRUE)
proxy5_dsP_WAPLS_RMSE <- sqrt(mean((df$DownscaledP - df$WA.PLS_prec)^2, na.rm = TRUE))


################# PLOTTING ################
all_coords <- data.frame(longitude=proxies$Longitude,latitude=proxies$Latitude, site = proxies$Site)
all_coords <- unique(all_coords)

coords <- SpatialPoints(cbind(all_coords$longitude, all_coords$latitude)) 

all_coords$ID <- c(4, 2, 1, 3, 5) # Order they should appear on map left to right


pastclim:::download_worldclim(var = "tavg", res = 10) # WorldClim data - Jan 
tavg_hres <- pastclim:::load_worldclim("tavg", 10)
tavg_hres <- terra::crop(tavg_hres, sea_ext)
tavg_hres_rast <- tavg_hres[[1]]


tavg_bio01 <- region_series(time_bp = 0, 
                            bio_variables = "bio01", 
                            dataset = 'Beyer2020', 
                           ext = sea_ext) # Raw model data

tdown_bio01 <- region_series(time_bp = 0, 
                            bio_variables = "bio01", 
                            dataset = 'custom', 
                            path_to_nc = "bio01_downscaled.nc") # Downscaled data

# data processing
library(ggplot2)
library(patchwork)
# spatial
library(sp)
library(raster)
library(rasterVis)
library(rgdal)
library(viridis)

colr <- colorRampPalette(rev(viridis(20)))



l1 <- levelplot(tavg_hres_rast, 
                margin=FALSE,                       
                colorkey=list(space='bottom',                   
                labels=list(at= seq(-50,40, 5), font = 4),
                axis.line=list(col='black')),    
                par.settings=list(axis.line=list(col='transparent'),
                panel.background = list(col= "lightgray")),
                scales=list(draw=FALSE),
                col.regions=colr, 
                at=c(seq(-45, 30, 2)),
                main = "WorldClim 1.2 January") + 
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 1.5)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)+5, c("4", "2", "1", "3", "5")))

l1

# Proxy 2
proxy2_ext <-  terra::ext(-91, -86, 36, 41) # 
proxy2_raw_plot <- crop(tavg_bio01$bio01, proxy2_ext)
proxy2_ds_plot <- crop(tdown_bio01$bio01, proxy2_ext)

l2 <- levelplot(proxy2_raw_plot, 
                margin=FALSE,                       
                colorkey=FALSE,  
                par.settings=list(axis.line=list(col='transparent'), 
                                  panel.background = list(col= "lightgray")),
                scales=list(draw=FALSE),            
                col.regions=colr, 
                at=c(seq(-45, 30, 5)),
                main = "Pittsburg Basin")+
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)-0.4, c("4", "2", "1", "3", "5"), cex = 2))


l3 <- levelplot(proxy2_ds_plot, 
                margin=FALSE,                       
                colorkey= FALSE, 
                par.settings=list(axis.line=list(col='transparent'),
                                  panel.background = list(col= "lightgray")),
                scales=list(draw=FALSE),            
                col.regions=colr, 
                at=c(seq(-45, 30, 1))) + 
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)-0.4, c("4", "2", "1", "3", "5"), cex = 2))


# Proxy 3
proxy3_ext <-  terra::ext(13, 17, 38, 42) # 
proxy3_raw_plot <- crop(tavg_bio01$bio01, proxy3_ext)
proxy3_ds_plot <- crop(tdown_bio01$bio01, proxy3_ext)

l4 <- levelplot(proxy3_raw_plot, 
                margin=FALSE,                       
                colorkey=FALSE,  
                par.settings=list(axis.line=list(col='transparent'), 
                                  panel.background = list(col= "lightgray")),
                scales=list(draw=FALSE),            
                col.regions=colr, 
                at=c(seq(-45, 30, 5)),
                main = "Lago Grande di Monticchio")+
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)-0.4, c("4", "2", "1", "3", "5"), cex = 2))


l5 <- levelplot(proxy3_ds_plot, 
                margin=FALSE,                       
                colorkey= FALSE, 
                par.settings=list(axis.line=list(col='transparent'),
                                  panel.background = list(col= "lightgray")),
                scales=list(draw=FALSE),            
                col.regions=colr, 
                at=c(seq(-45, 30, 1))) + 
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)-0.4, c("4", "2", "1", "3", "5"), cex = 2))


# Proxy 4

proxy4_ext <-  terra::ext(99, 103, 36, 40) # 
proxy4_raw_plot <- crop(tavg_bio01$bio01, proxy4_ext)
proxy4_ds_plot <- crop(tdown_bio01$bio01, proxy4_ext)

l6 <- levelplot(proxy4_raw_plot, 
                margin=FALSE,                       
                colorkey=FALSE,  
                par.settings=list(axis.line=list(col='transparent')),
                scales=list(draw=FALSE),  
                col.regions=colr, 
                at=c(seq(-45, 30, 2)),
                main = "Luanhaizi LH2")+
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)+0.4, c("4", "2", "1", "3", "5"), cex = 2))


l7 <- levelplot(proxy4_ds_plot, 
                margin=FALSE,                       
                colorkey=FALSE,   
                par.settings=list(axis.line=list(col='transparent')),
                scales=list(draw=FALSE),            
                col.regions=colr, 
                at=c(seq(-45, 30, 2)))+
  latticeExtra:::layer(sp.points(coords, col = "black", pch = 20, cex = 2)) +
  latticeExtra:::layer(sp.text(cbind(all_coords$longitude, all_coords$latitude)+0.4, c("4", "2", "1", "3", "5"), cex = 2))

# DOESN'T WORK

require(gridExtra) # also loads grid
require(lattice)
gridExtra:::grid.arrange(l2, l4, l6, l3, l5, l7, ncol =3, nrow = 2)


