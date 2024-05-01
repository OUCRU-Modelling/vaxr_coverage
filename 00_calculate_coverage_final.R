
# read data ---------------------------------------------------------------

library(tidyverse)
library(arrow)

## cleaned personal info
data_path_children <- file.path("/.", "cluster_data", "vrdata", "cleaned", "personal_info.parquet") 
children <- open_dataset(data_path_children)

## cleaned vacname data
vacname_data_path <- file.path("/.", "cluster_data", "vrdata", "cleaned", "vacname") 
vacname_data <- open_dataset(vacname_data_path)
children_province_reg <- vacname_data %>%
  collect() %>%
  group_by(pid) %>%
  arrange(vacdate) %>%
  summarise(province_reg = province_reg[1]) %>%
  ungroup() 

## cleaned pathogen data
data_path_pathogen <- file.path("/.", "cluster_data", "vrdata", "cleaned", "pathogen") 
pathogen_data <- open_dataset(data_path_pathogen)
pathogen <- as.vector(unique(pathogen_data %>% select(pathogen) %>% collect()))

## coverage by province & age group
children_final <- merge(children, children_province_reg, by = "pid") %>%
  mutate(yob = year(dob),
         province_doc = province,
         province_vac = province_reg) %>%
  filter(!province_doc %in% c("Tỉnh tập huấn")) %>%
  select(pid, sex, dob, yob, province_doc, province_vac)

children_sum_doc <- children_final %>%
  group_by(yob, province_doc) %>%
  summarise(n_doc = n()) %>%
  ungroup() %>%
  mutate(cohort_doc = 1:n())

children_sum_vac <- children_final %>%
  group_by(yob, province_vac) %>%
  summarise(n_vac = n()) %>%
  ungroup() %>%
  mutate(cohort_vac = 1:n())

children_sum_combine <- rbind(
  children_sum_doc %>%
    rename(province = province_doc,
           n = n_doc,
           cohort = cohort_doc) %>%
    mutate(cohort_type = "cohort_doc")
)

children_final_add <- merge(merge(children_final, children_sum_doc, by = c("yob", "province_doc"), all.x = TRUE),
                            children_sum_vac, by = c("yob", "province_vac"), all.x = TRUE)
save(children_final, children_sum_doc, children_sum_vac, children_sum_combine, children_final_add, file = file.path("save", "children.Rdata"))

get_coverage_agegroup <- function(pathogen) {
  data_path_tmp <- file.path(data_path_pathogen, paste0("pathogen=", pathogen))
  tmp <- open_dataset(data_path_tmp)
  
  children_path <- merge(tmp, children_final_add, by = "pid", all.x = TRUE)
  
  children_path_sum <- children_path %>%
    mutate(vacage = (vacdate - dob)/dyears(1),
           vacage_group = ifelse(vacage <= 1, "(0-1y]",
                                 ifelse(vacage <= 2, "(1-2y]",
                                        ifelse(vacage <= 3, "(2-3]", "(3+)")))) %>%
    group_by(pid, vacage_group) %>%
    arrange(vacage) %>%
    summarise(cohort_doc = cohort_doc[1],
              cohort_vac = cohort_vac[1],
              n_shot = n(),
              first_shot_age = vacage[1]) %>%
    ungroup()
  
  children_path_sum_all <- children_path_sum %>%
    group_by(pid) %>%
    summarise(cohort_doc = cohort_doc[1],
              cohort_vac = cohort_vac[1],
              n_shot = sum(n_shot, na.rm = TRUE),
              first_shot_age = min(first_shot_age, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(vacage_group = "all") %>%
    select(pid, vacage_group, cohort_doc, cohort_vac, n_shot, first_shot_age)
  
  children_path_sum_combine <- rbind(children_path_sum, children_path_sum_all) %>%
    mutate(vacage_group = factor(vacage_group, levels = c("(0-1y]", "(1-2y]", "(2-3]", "(3+)", "all")))
  
  children_path_sum_final_doc <- children_path_sum_combine %>%
    group_by(cohort_doc, vacage_group) %>%
    summarise(anyshot   = sum(n_shot > 0),
              shots_1     = sum(n_shot == 1),
              shots_2     = sum(n_shot == 2),
              shots_3     = sum(n_shot == 3),
              shots_4     = sum(n_shot == 4),
              first_shot_age_min = min(first_shot_age),
              first_shot_age_mean = mean(first_shot_age),
              first_shot_age_median = median(first_shot_age),
              first_shot_age_max = max(first_shot_age)) %>%
    ungroup()
  
  children_path_sum_final_vac <- children_path_sum_combine %>%
    group_by(cohort_vac, vacage_group) %>%
    summarise(anyshot   = sum(n_shot > 0),
              shots_1     = sum(n_shot == 1),
              shots_2     = sum(n_shot == 2),
              shots_3     = sum(n_shot == 3),
              shots_4     = sum(n_shot == 4),
              first_shot_age_min = min(first_shot_age),
              first_shot_age_mean = mean(first_shot_age),
              first_shot_age_median = median(first_shot_age),
              first_shot_age_max = max(first_shot_age)) %>%
    ungroup()
  
  output_doc <- merge(children_sum_doc, children_path_sum_final_doc, by = "cohort_doc", all.x = TRUE) %>%
    rename(province = province_doc,
           n = n_doc,
           cohort = cohort_doc) %>%
    mutate(cohort_type = "cohort_doc")
  
  output_vac <- merge(children_sum_vac, children_path_sum_final_vac, by = "cohort_vac", all.x = TRUE) %>%
    rename(province = province_vac,
           n = n_vac,
           cohort = cohort_vac) %>%
    mutate(cohort_type = "cohort_vac")
  
  output <- rbind(output_doc, output_vac) %>%
    mutate(coverage = anyshot/n,
           pathogen = pathogen)
  
  save(output, file = file.path("save", paste0("coverage_", pathogen, ".Rdata")))
  return(output)
}

### get coverage 
coverage_province_vacage <- do.call(what = "rbind", 
                                    args = lapply(X = pathogen$pathogen, 
                                                  FUN = get_coverage_agegroup)) 

coverage_vacage <- coverage_province_vacage %>%
  group_by(yob, pathogen, vacage_group) %>%
  summarise(
    n = sum(n),
    anyshot = sum(anyshot, na.rm = TRUE),
    shots_1     = sum(shots_1, na.rm = TRUE),
    shots_2     = sum(shots_2, na.rm = TRUE),
    shots_3     = sum(shots_3, na.rm = TRUE),
    shots_4     = sum(shots_4, na.rm = TRUE),
    coverage = anyshot/n
  )

save(coverage_province_vacage, coverage_vacage, file = file.path("save", "coverage_sum.Rdata"))

## for dashboard
db_coverage_province_vacage <- coverage_province_vacage %>%
  filter((pathogen %in% c("bcg", "dip", "hepb", "hib", "mum", "nm", "mea", "jev", "rub", "spn")) & 
           (!is.na(vacage_group) & vacage_group %in% c("(0-1y]", "all"))) %>%
  mutate(pathogen = case_when(
    pathogen == "bcg" ~ "Tuberculosis",
    pathogen == "dip" ~ "Diptheria",
    pathogen == "hepb" ~ "Hepatitis B",
    pathogen == "hib" ~ "Haemophilus influenzae type b",
    pathogen == "mum" ~ "Mumps",
    pathogen == "nm" ~ "Meningococcus",
    pathogen == "mea" ~ "Measles",
    pathogen == "jev" ~ "Japanese Encephalitis",
    pathogen == "rub" ~ "Rubella",
    pathogen == "spn" ~ "Streptococcus pneumoniae",
  )) %>%
  mutate(province = ifelse(province == "Thành phố Hồ Chí Minh", "Hồ Chí Minh",
                           ifelse(province == "Hòa Bình", "Hoà Bình", province))) %>%
  select(pathogen, yob, province, n, vacage_group, anyshot, shots_1, shots_2, shots_3, shots_4, coverage)

db_coverage_vacage <- coverage_vacage %>%
  filter((pathogen %in% c("bcg", "dip", "hepb", "hib", "mum", "nm", "mea", "jev", "rub", "spn")) & 
           (!is.na(vacage_group) & vacage_group %in% c("(0-1y]", "all"))) %>%
  mutate(pathogen = case_when(
    pathogen == "bcg" ~ "Tuberculosis",
    pathogen == "dip" ~ "Diptheria",
    pathogen == "hepb" ~ "Hepatitis B",
    pathogen == "hib" ~ "Haemophilus influenzae type b",
    pathogen == "mum" ~ "Mumps",
    pathogen == "nm" ~ "Meningococcus",
    pathogen == "mea" ~ "Measles",
    pathogen == "jev" ~ "Japanese Encephalitis",
    pathogen == "rub" ~ "Rubella",
    pathogen == "spn" ~ "Streptococcus pneumoniae",
  )) %>%
  select(pathogen, yob, n, vacage_group, anyshot, shots_1, shots_2, shots_3, shots_4, coverage)

save(db_coverage_province_vacage, db_coverage_vacage, file = file.path("save", "db_coverage_sum.Rdata"))

write_csv(db_coverage_province_vacage, file = file.path("save", "db_coverage_province_vacage.csv"))
write_csv(db_coverage_vacage, file = file.path("save", "db_coverage_vacage.csv"))

## plot
library(ggplot2)
ggplot(data = coverage_all_sum, aes(x = yob)) +
  geom_line(aes(y = all_coverage)) +
  geom_line(data = coverage_all_sum2, aes(y = all_coverage), color = "red") +
  facet_wrap(~ pathogen, scales = "free_y") +
  theme_bw()
ggsave(filename = file.path("figures", "coverage_all_sum.png"), width = 15, height = 10)


# map ---------------------------------------------------------------------

#install.packages("terra")
#remotes::install_github("r-spatial/sf")
#install.packages("maptools", repos="http://R-Forge.R-project.org")
#remotes::install_github('oswaldosantos/ggsn')

library(ggplot2)
library(sf)
library(ggsn) 

vn_tinh <- st_read(dsn = file.path("save", "gadm41_VNM.gpkg"), layer = "ADM_ADM_1")
vn_qh <- st_read(dsn = file.path(map_path, "gadm41_VNM.gpkg"), layer = "ADM_ADM_2")
vn_px <- st_read(dsn = file.path(map_path, "gadm41_VNM.gpkg"), layer = "ADM_ADM_3")

## match province names
province_name <- data.frame(province = unique(coverage_all2$province),
                            NAME_1 = vn_tinh$NAME_1[match(unique(coverage_all2$province), vn_tinh$NAME_1)])
province_name$NAME_1[province_name$province == "Thành phố Hồ Chí Minh"] <- "Hồ Chí Minh"
province_name$NAME_1[province_name$province == "Hòa Bình"] <- "Hoà Bình"

coverage_all3 <- merge(merge(coverage_all2, province_name, by = "province", all.x = TRUE), vn_tinh, by = "NAME_1", all.x = TRUE)
save(coverage_all2, vn_tinh, province_name, coverage_all3, file = file.path("save", "coverage_all3.Rdata"))  

## xem dữ liệu trong vn_tinh
str(vn_tinh)
vn_tinh$NAME_1

## chọn bản đồ TPHCM cấp tỉnh
hcm_shp <- vn_tinh %>%
  filter(NAME_1 == "Hồ Chí Minh")

## chọn bản đồ TPHCM cấp quận/huyện
hcm_shp_qh <- vn_qh %>%
  filter(NAME_1 == "Hồ Chí Minh")

## chọn bản đồ TPHCM cấp phường/xã
hcm_shp_px <- vn_px %>%
  filter(NAME_1 == "Hồ Chí Minh")

## vẽ thử
hcm_shp %>%
  ggplot() + geom_sf()

hcm_shp_qh %>%
  ggplot() + geom_sf()

hcm_shp_px %>%
  ggplot() + geom_sf()

## lấy thông tin mã bản đồ
hcm_map_id <- st_drop_geometry(hcm_shp_px) %>%
  transmute(
    qh = NAME_2,
    px = NAME_3,
    GID_2 = GID_2,
    GID_3 = GID_3
  )
View(hcm_map_id)

tmp <- merge(hcm_map_id, 
             ma_px_moi %>%
               dplyr::select(qh, px, ma_px_2022),
             by = c("qh", "px"),
             all = TRUE)
write_xlsx(x = tmp, path = file.path(map_path, "map_id.xlsx"))

## tạo bản đồ cấp phường xã và quận huyện mới sau khi sáp nhập
map_id_final <- read_excel(path = file.path(map_path, "map_id_final.xlsx"))
map_id_final_qh <- map_id_final %>%
  dplyr::select(qh, GID_2)

### cấp quận huyện
hcm_shp_qh2 <- left_join(hcm_shp_qh, map_id_final_qh, by = c("GID_2" = "GID_2")) %>%
  group_by(qh) %>%
  summarise(geometry = st_union(geom)) %>%
  ungroup()

###- kiểm tra
hcm_shp_qh2 %>%
  ggplot() +
  geom_sf() +
  geom_sf_label(aes(label = qh))

### cấp phường xã
hcm_shp_px2 <- left_join(hcm_shp_px, map_id_final, by = c("GID_2" = "GID_2", "GID_3" = "GID_3")) %>%
  group_by(qh, px, ma_px_2022) %>%
  summarise(geometry = st_union(geom)) %>%
  ungroup()

###- kiểm tra
hcm_shp_px2 %>%
  ggplot() +
  geom_sf() +
  geom_sf_label(aes(label = ma_px_2022))

## lưu dữ liệu
save(hcm_shp, hcm_shp_qh2, hcm_shp_px2, file = file.path(rdata_path, "hcm_map.Rdata"))

qh_tiemchung <- left_join(hcm_shp_qh2, muitiem_tile_qh, by = c("qh" = "qh"))
px_tiemchung <- left_join(hcm_shp_px2, muitiem_tile_px, by = c("qh" = "qh", "px" = "px", "ma_px_2022" = "ma_px_2022"))

## chọn cách chia màu (https://colorbrewer2.org)
range(qh_tiemchung$p_tcdd)
mypal <- function(x) c("#e5f5e0", "#a1d99b", "#31a354")
mybreak <- c(85, 90, 95, 100)


# new ---------------------------------------------------------------------

file.copy(from = file.path("data"), to = file.path("docs"), 
          overwrite = TRUE, recursive = TRUE, 
          copy.mode = TRUE)
file.copy(from = file.path("dashboard_files"), to = file.path("docs"), 
          overwrite = TRUE, recursive = TRUE, 
          copy.mode = TRUE)
