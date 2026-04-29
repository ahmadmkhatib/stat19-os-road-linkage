


# create 500m buffer around dissolved CAZ
caz_buffer_500 <- st_buffer(caz_dissolved, 500)

# roads inside buffer
roads_in_buffer <- st_intersects(road_attributes, caz_buffer_500, sparse = FALSE) %>%
  as_tibble() %>%
  mutate(identifier = road_attributes$identifier) %>%
  filter(value == TRUE) %>%
  select(identifier)

# flag buffer roads
road_buffer_flag <- road_attributes %>%
  mutate(
    in_buffer = identifier %in% roads_in_buffer$identifier
  ) %>%
  st_drop_geometry() %>%
  select(identifier, in_buffer)


#Explicit Spillover Ring Variable


road_attributes <- road_attributes %>%
  mutate(
    inside_caz = identifier %in% road_caz_roadlevel$identifier
  )

buffer_500  <- st_buffer(caz_dissolved, 500)
buffer_1000 <- st_buffer(caz_dissolved, 1000)

road_attributes <- road_attributes %>%
  mutate(
    ring_500  = st_intersects(geometry, buffer_500, sparse = FALSE)[,1],
    ring_1000 = st_intersects(geometry, buffer_1000, sparse = FALSE)[,1]
  )


#  join to panel and:

analysis_data <- road_panel_complete %>%
  left_join(road_buffer_flag, by = "identifier") %>%
  filter(in_buffer == FALSE | ever_treated == 1)

