coords <- jsonlite::read_json("analysis/data/raw/ex_coords.txt", simplifyVector = TRUE)

# load packages
library(ggplot2)
library(rworldmap)
library(viridisLite)
library(gridExtra)

# -----------------------------------
## STEP 1: MAKE A BASIC MAP
# -----------------------------------


plot_df <- cbind(coords, seq_len(nrow(coords)))
names(plot_df) <- c("lat", "long", "measure")

# get into long format for facetted plot
plot_df_long <- tidyr::gather(data = plot_df, key = component, value = value, 3, factor_key = TRUE)

# -----------------------------------
# PLOTTING

# define plotting parameters
col_country <- grey(0.95)
col_country_border <- grey(0.5)
size_country_border <- 0.5
col_sea <- grey(1.0)
resolution <- "high"
col_limits <- range(plot_df$measure)
point_size <- 4
stroke_col <- grey(0.5)
stroke_size <- 0.25
col_vec <- rev(colorRampPalette(bobfunctions2::col_tim())(nrow(plot_df)))

# load country shapefiles
world_map <- getMap(resolution = resolution)

# basic map plot
plot_base <- ggplot() + theme_bw() + theme(panel.background = element_rect(fill = col_sea),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           strip.background = element_blank(),
                                           strip.text = element_text(angle = 0, hjust = 0, size = 12))

# add country borders
plot_base <- plot_base + geom_polygon(aes(long, lat, group = group),
                                      size = size_country_border, color = col_country_border,
                                      fill = col_country, data = world_map)


# make separate inset plots
locations_gg <- plot_base + coord_cartesian(xlim = c(28.5, 31), ylim = c(-3,-1)) +
  geom_point(aes(x = long, y = lat, fill = as.factor(value)), color = stroke_col, shape = 21,
             stroke = stroke_size, size = point_size, data = plot_df_long) +
  scale_fill_manual(values = col_vec, name = "Location") +
  xlab("longitude") + ylab("latitude") +
  ggtitle("Locations of Sample Sites")

# -----------------------------------
## STEP 2: Simulate some data
# -----------------------------------

base_df <- plot_df[,1:2]
base_df$R561H <- 0

# Simulate our radiation
# break and distribute counts
centre <- 8:20
sw <- 1:2
n <- 3:6
ne <- 7
s <- 21:26
se <- 27

base_df$R561H[centre] <- runif(length(centre), min = 0.1, max = 0.3)
base_df$R561H[sw] <- runif(length(sw), min = 0, max = 0.02)
base_df$R561H[n] <- runif(length(n), min = 0, max = 0.1)
base_df$R561H[ne] <- runif(length(ne), min = 0, max = 0.05)
base_df$R561H[s] <- runif(length(s), min = 0.1, max = 0.15)
base_df$R561H[se] <- runif(length(se), min = 0, max = 0.04)

# simulate our N
base_df$N1 <- 100
base_df$N2 <- 25
base_df$N2[c(centre, s)] <- 200
base_df$N3 <- 25
base_df$N3[c(sw)] <- 150
base_df$N3[c(n, ne, 8:11)] <- 25
base_df$N3[c(12:20)] <- 150
base_df$N3[c(21:22)] <- 200
base_df$N3[c(23:26)] <- 150

# base_df$N <- extraDistr::rbbinom(seq_along(base_df$lat), size = 100, 1, 0.1)

base_df$R561H_N1 <- rbinom(length(base_df$N1), base_df$N1, base_df$R561H)
base_df$R561H_N2 <- rbinom(length(base_df$N2), base_df$N2, base_df$R561H)
base_df$R561H_N3 <- rbinom(length(base_df$N3), base_df$N3, base_df$R561H)

# logit transform
eps <- 0.5
base_df$R561H_N1_logit <- log((base_df$R561H_N1 + eps)/(base_df$N1 - base_df$R561H_N1 + eps))
base_df$R561H_N2_logit <- log((base_df$R561H_N2 + eps)/(base_df$N2 - base_df$R561H_N2 + eps))
base_df$R561H_N3_logit <- log((base_df$R561H_N3 + eps)/(base_df$N3 - base_df$R561H_N3 + eps))

# fit prevmap
fit <- list()
fit$R561H_N1_fit <- PrevMap::linear.model.MLE(formula = R561H_N1_logit ~ 1, coords = ~long + lat,
                                          data = base_df, start.cov.pars = c(3,3), kappa = 1)
fit$R561H_N2_fit <- PrevMap::linear.model.MLE(formula = R561H_N2_logit ~ 1, coords = ~long + lat,
                                           data = base_df, start.cov.pars = c(3,3), kappa = 1)
fit$R561H_N3_fit <- PrevMap::linear.model.MLE(formula = R561H_N3_logit ~ 1, coords = ~long + lat,
                                           data = base_df, start.cov.pars = c(3,3), kappa = 1)

# create grid of points to map
poly <- as.matrix(expand.grid(c(28.5, 30, 31.5), c(-3,-2, -1)))[-5,]
poly <- poly[order(poly[,1]),]
poly <- poly[c(1:3,5,8:6,4),]
grid_pred <- splancs::gridpts(poly, xs = 0.05, ys = 0.05)
colnames(grid_pred) <- c("long","lat")

# make list of model predictions from fits
pred <- list()
for (i in 1:length(fit)) {

  # make model predictions
  pred_raw <- PrevMap::spatial.pred.linear.MLE(fit[[i]], grid_pred, scale.predictions = "prevalence",
                                      n.sim.prev = 1e3, standard.errors = TRUE)

  # make predictions into raster then dataframe
  rast <- raster::rasterFromXYZ(cbind(pred_raw$grid.pred, pred_raw$prevalence$prediction))
  df_rast <- as.data.frame(raster::rasterToPoints(rast))

  rast_low <- raster::rasterFromXYZ(cbind(pred_raw$grid.pred, pred_raw$prevalence$quantiles[,1]))
  rast_high <- raster::rasterFromXYZ(cbind(pred_raw$grid.pred, pred_raw$prevalence$quantiles[,2]))
  df_rast_low <- as.data.frame(raster::rasterToPoints(rast_low))
  df_rast_high <- as.data.frame(raster::rasterToPoints(rast_high))

  df_rast <- cbind(df_rast, df_rast_low[,3], df_rast_high[,3])
  names(df_rast)[4:5] <- c("low", "high")

  # make smooth raster (for use when producing contour lines)
  rast_smooth <- raster::disaggregate(raster::aggregate(rast, fact = 3), fact = 3, method = "bilinear")
  rast_smooth_low <- raster::disaggregate(raster::aggregate(rast_low, fact = 3), fact = 3, method = "bilinear")
  rast_smooth_high <- raster::disaggregate(raster::aggregate(rast_high, fact = 3), fact = 3, method = "bilinear")
  df_rast_smooth <- as.data.frame(raster::rasterToPoints(rast_smooth))
  df_rast_smooth_low <- as.data.frame(raster::rasterToPoints(rast_smooth_low))
  df_rast_smooth_high <- as.data.frame(raster::rasterToPoints(rast_smooth_high))

  df_rast_smooth <- cbind(df_rast_smooth, df_rast_smooth_low[,3], df_rast_smooth_high[,3])
  names(df_rast_smooth)[4:5] <- c("low", "high")

  # store in pred list
  pred[[i]] <- list(df_rast = df_rast,
                    df_rast_smooth = df_rast_smooth)
}

# save to file
saveRDS(pred, "analysis/data/derived/ex_prevmap.rds")

# -----------------------------------
## STEP 3: Plot Our predictions
# -----------------------------------

# define plotting parameters
col_country <- grey(0.95)
col_country_border <- grey(0.5)
size_country_border <- 0.5
#col_sea <- grey(1.0)
shape_resolution <- "high"
#col_limits <- c(-4,4)
#point_size <- 2.5
#stroke_col <- grey(0.5)
#stroke_size <- 0.25
col_vec <- viridisLite::plasma(100)

# read in predictions
pred <- readRDS("analysis/data/derived/ex_prevmap.rds")

# load country shapefiles, leaving gap for DRC
world_map <- getMap(resolution = shape_resolution)
sub_world_map <- subset(world_map, ISO3 != "RWA")

# create list of plot objects
plot_list_r561h_sample_plots <- list()

for(i in seq_along(pred)) {

plot_list_r561h <- list()
plot_list_r561h[[1]] <- plot_base +
  geom_raster(aes(x = x, y = y, fill = layer), data = pred[[i]]$df_rast) +
  geom_contour(aes(x = x, y = y, z = layer), size = 0.5, color = "black", data = pred[[i]]$df_rast_smooth) +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = sub_world_map) +
  coord_cartesian(xlim = c(28.5, 31), ylim = c(-3,-1)) +
  scale_fill_gradientn(colours = col_vec, name = "prevalence", limits = c(0,0.5)) +
  xlab("longitude") + ylab("latitude") +
  ggtitle(expression("R561H prevalence"))

plot_list_r561h[[2]] <- plot_base +
  geom_raster(aes(x = x, y = y, fill = low), data = pred[[i]]$df_rast) +
  geom_contour(aes(x = x, y = y, z = low), size = 0.5, color = "black", data = pred[[i]]$df_rast_smooth) +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = sub_world_map) +
  coord_cartesian(xlim = c(28.5, 31), ylim = c(-3,-1)) +
  scale_fill_gradientn(colours = col_vec, name = "prevalence", limits = c(0,0.5)) +
  xlab("longitude") + ylab("latitude") +
  ggtitle(expression("R561H prevalence 2.5% CI"))

plot_list_r561h[[3]] <- plot_base +
  geom_raster(aes(x = x, y = y, fill = high), data = pred[[i]]$df_rast) +
  geom_contour(aes(x = x, y = y, z = high), size = 0.5, color = "black", data = pred[[i]]$df_rast_smooth) +
  geom_polygon(aes(long, lat, group = group),
               size = size_country_border, color = col_country_border,
               fill = col_country, data = sub_world_map) +
  coord_cartesian(xlim = c(28.5, 31), ylim = c(-3,-1)) +
  scale_fill_gradientn(colours = col_vec, name = "prevalence", limits = c(0,0.5)) +
  xlab("longitude") + ylab("latitude") +
  ggtitle(expression("R561H prevalence 97.5% CI"))


# finalise plots
prev_map_r561h <- cowplot::plot_grid(plot_list_r561h[[2]], plot_list_r561h[[1]], plot_list_r561h[[3]], ncol = 3)
plot_list_r561h_sample_plots[[i]] <- prev_map_r561h

}

# -----------------------------------
## STEP 4: Plot The Observed Data
# -----------------------------------

# make separate inset plots
raw_data_gg <- plot_base + coord_cartesian(xlim = c(28.5, 31), ylim = c(-3,-1)) +
  geom_point(aes(x = long, y = lat, fill = R561H/N), color = stroke_col, shape = 21,
             stroke = stroke_size, size = point_size, data = base_df) +
  scale_fill_gradientn(colours = col_vec, name = "prevalence", limits = c(0,0.3)) +
  xlab("longitude") + ylab("latitude") +
  ggtitle("Observed 561H Prevalence at Sample Sites")


tr1 <- cowplot::plot_grid(NULL, locations_gg, NULL, raw_data_gg, NULL, rel_widths = c(1/6-0.025,1/3,0.05, 1/3,1/6-0.025), nrow = 1)

comb1 <- cowplot::plot_grid(tr1, plot_list_r561h_sample_plots[[1]], ncol = 1)
comb2 <- cowplot::plot_grid(tr1, plot_list_r561h_sample_plots[[2]], ncol = 1)
comb3 <- cowplot::plot_grid(tr1, plot_list_r561h_sample_plots[[3]], ncol = 1)

save_figs("dummy1", comb1, 15, 8)
save_figs("dummy2", comb2, 15, 8)
save_figs("dummy3", comb3, 15, 8)
