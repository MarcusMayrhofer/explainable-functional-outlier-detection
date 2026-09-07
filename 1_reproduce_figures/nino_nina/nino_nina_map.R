###############################################################################
# Map of the four Niño regions in the Pacific Ocean (figure: plots/map.pdf) (created using claude code)
#
# Region definitions (source: https://www.cpc.ncep.noaa.gov/data/indices/, ERSSTv5):
#   Niño 1+2  (0-10°S, 90°W-80°W)
#   Niño 3    (5°N-5°S, 150°W-90°W)
#   Niño 4    (5°N-5°S, 160°E-150°W)
#   Niño 3.4  (5°N-5°S, 170°W-120°W)
###############################################################################

library(tidyverse)
library(sf)
library(rnaturalearth)
require(this.path)
setwd(this.path::this.dir())

# World polygons
worldMap <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid()

# Projection centred on the Pacific (longitude 133°E)
target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=133")

# Thin polygon along lon = 180 - 133 used to slice the world
# Polygons so that countries straddling the projection seam are not drawn across it
offset <- 180 - 133
polygon <- st_polygon(x = list(rbind(
  c(-0.0001 - offset, 90),
  c(0 - offset, 90),
  c(0 - offset, -90),
  c(-0.0001 - offset, -90),
  c(-0.0001 - offset, 90)
))) %>%
  st_sfc() %>%
  st_set_crs(4326)

# Subtract the seam polygon and reproject to the Pacific-centred map
world2 <- worldMap %>% st_difference(polygon)
world3 <- world2 %>% st_transform(crs = target_crs)

# Crop to the Pacific bounding box used for the figure
world4 <- st_crop(
  x = world3,
  y = st_as_sfc(
    st_bbox(c(xmin = -60, xmax = 120, ymin = -40, ymax = 40), crs = 4326)
  ) %>% st_transform(target_crs)
)

# Bounding boxes of the four Niño regions.
rectangles <- data.frame(
  xmin = c(-90, -150, -150, -170),
  xmax = c(-80, -90, 160, -120),
  ymin = c(-10, -5, -5, -5),
  ymax = c(0, 5, 5, 5),
  fill = c("Niño 1+2", "Niño 3", "Niño 4", "Niño 3.4"),
  linetype = c("A", "A", "A", "B"),
  label = c("Niño 1+2", "Niño 3", "Niño 4", "Niño 3.4")
)

# Convert each row to an sf polygon and reproject to the map
rectangles_sf <- rectangles %>%
  mutate(geometry = pmap(
    list(xmin, xmax, ymin, ymax), 
    ~ st_polygon(list(rbind(
      c(..1, ..3),
      c(..2, ..3),
      c(..2, ..4),
      c(..1, ..4),
      c(..1, ..3)
    )))
  )) %>%
  st_as_sf(crs = 4326) %>%
  st_transform(target_crs)

# Label positions: centroid of each rectangle, shifted vertically so labels do not overlap
label_positions <- rectangles_sf %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  setNames(c("x", "y")) %>%
  mutate("labels" = rectangles_sf$label, linetype = rectangles_sf$linetype)
label_positions$y <- label_positions$y + c(-900000, 900000, 900000, -900000)

rectangles_sf <- rectangles_sf %>%
  mutate(x = label_positions$x, y = label_positions$y)


# Final map: world coastline in grey, Niño regions as filled rectangles with text labels
pacific_plot <- ggplot(data = world4) +
  geom_sf(fill = "grey") +
  geom_sf(data = rectangles_sf, aes(geometry = geometry, fill = fill, 
                                    linetype = linetype, color = linetype), linewidth = 0.7, alpha = 0.4) +
  geom_text(data = rectangles_sf, aes(x = x, y = y, label = label, color = linetype), size = 5) +
  coord_sf(crs = target_crs) +
  scale_color_manual(values = c("#222222", "#bf0000"), guide = "none") + 
  scale_fill_manual(values = c("darkgreen", "darkorange", "transparent", "darkblue"), guide = "none") + 
  scale_x_continuous(expand = c(0,0), breaks = seq(from = -180, to = 180, by = 20)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_linetype_manual(values = c(1,2), guide = "none") + 
  theme_bw() + 
  labs(x = element_blank(), y = element_blank())
pacific_plot
ggsave(plot = pacific_plot, filename = "plots/map.pdf", width = 9, height = 4)
