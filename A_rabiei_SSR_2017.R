devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "util.R")
# install.packages("Rcpp")
# library(Rcpp)

# devtools::install_github("tidyverse/readxl")
# devtools::install_github("tidyverse/dplyr")
# devtools::install_github("tidyverse/tidyverse")
# devtools::install_github("tidyverse/tidyverse")
pacman::p_load(char = c("gtools","poppr",  "tidyverse", "bindr", "ape", "DataExplorer", 
                        "readxl", "scales"))
# .fixdevtools()

#### Load raw data from AGRF ####
AGRF_sample_file <- recent_file("./data", "2017 POP GEN_DNA Extraction.*.xlsx")
# AGRF_batches <- list.dirs("..", recursive = FALSE) %>% .[grepl("QAGRF", .)] %>% str_extract(., "QAGRF.+")

# Load sample processing (extraction and sequencing) information
# AGRF_samples <- excel_sheets(AGRF_sample_file) %>% .[grepl("ABMS", .)] %>%
#   map_dfr(~{read_excel(AGRF_sample_file, sheet = ., cell_cols("A:N")) %>% 
#       mutate(AGRF_Batch=sub("_samples", "", .x), AGRF_Tube_Id=as.character(AGRF_Tube_Id)) %>% select(AGRF_Tube_Id, Sample_name, AGRF_Batch)} )
# Load isolate information
AGRF_samples <- read_excel(AGRF_sample_file, sheet = "Sample_information") %>% 
  mutate(AGRF_Tube_Id=as.character(AGRF_Tube_Id))
  # filter(!grepl("remove|ignore", Comment, ignore.case = TRUE))

AGRF_files <- list.files("./data", "ABMS.+.xlsx", include.dirs = FALSE, 
                         full.names = TRUE, recursive = TRUE)

#### identify heterozygotes ####
# Check for non matching alleles (genotyping error/contamination/heterozugotes?)
heterozyg_samples <- AGRF_files %>% map_dfr(~read_excel(., sheet = "Genotypes") %>% 
                                  setNames(., gsub(" ", "", colnames(.))) %>% 
                                  filter_at(vars(Allele1:Allele2), all_vars(!grepl("Failed", .))) %>% 
                                  mutate_all(as.character) %>% 
                                  mutate(AGRF_Batch=str_extract(.x, "ABMS-[0-9]+"))) %>% 
  filter(!is.na(Allele1), !is.na(Allele2), Allele1!=Allele2) %>% 
  dplyr::select(c("SampleName", "AGRF_Batch")) %>% distinct()
  
  
  
  
# Save shterozygote amples to file
  AGRF_files %>% map_dfr(~read_excel(., sheet = "Genotypes") %>% 
                         setNames(., gsub(" ", "", colnames(.))) %>% 
                         filter_at(vars(Allele1:Allele2), all_vars(!grepl("Failed", .))) %>% 
    mutate_all(as.character) %>% mutate(AGRF_Batch=str_extract(.x, "ABMS-[0-9]+"))) %>% 
    inner_join(heterozyg_samples) %>% 
  dplyr::select(c("SampleName", "Marker", "Allele1", "Allele2", "AGRF_Batch")) %>% 
  gather(variable, value, Allele1:Allele2) %>% 
  unite(temp, Marker, variable) %>% spread(temp, value) %>% 
  inner_join(AGRF_samples %>% select(-starts_with("Extraction")) %>% 
               filter(!grepl("remove|ignore", Comment, ignore.case = TRUE)) %>% 
               select(-Comment), ., 
             by=c("AGRF_Tube_Id"="SampleName", "AGRF_Batch")) %>% # mutate(SampleName=Sample_Name) %>%
  select(-starts_with("AGRF")) %>% write_csv("./output/A_rabiei_SSR_2017_heterozygous_samples.csv")


SSR_data_2017_wide <- data.frame()
for (f in AGRF_files) {
  # f=AGRF_files[1]
  SSR_data_2017_wide <- rbind(SSR_data_2017_wide, read_excel(f, sheet = "Genotypes") %>% 
    setNames(., gsub(" ", "", colnames(.))) %>% mutate_all(as.character) %>% 
      mutate_at(vars(Allele1:Allele2), list(~if_else(grepl("Failed", .), NA_character_, .))) %>%
      # mutate(SampleName=as.character(SampleName)) %>% # filter
    dplyr::filter(!is.na(Allele1), Allele1==Allele2) %>% dplyr::select(1:3) %>% 
      # gather(variable, value, Allele1:Allele2) %>% 
      # unite(temp, Marker, variable) 
      spread(Marker, Allele1) %>% 
    inner_join(AGRF_samples %>% select(-starts_with("Extraction")) %>% 
                 filter(AGRF_Batch==str_extract(f, "ABMS-[0-9]+"), 
                        !grepl("remove|ignore", Comment, ignore.case = TRUE)) %>% 
                 select(-Comment), ., 
                 by=c("AGRF_Tube_Id"="SampleName")) %>% # mutate(SampleName=Sample_Name) %>%
    select(-starts_with("AGRF")))
}

# Number of valid genotyped samples
# SSR_data_2017_wide %>% filter(rowSums(is.na(.))<2) %>% nrow()
# remove missing data
# SSR_data_2017_wide <- SSR_data_2017_wide %>% # rename(Sample_Name=SampleName) %>% 
#   filter(rowSums(is.na(.))<2)


# Change NAs to 0s
SSR_data_2017_wide[is.na(SSR_data_2017_wide)] <- "0"

# Append to 4-year database file
# Read in sample details - Sam to add the following fields: Regions	isolate_name	location	State	Host	Year

SSR_data_2017_wide %>% write_csv("./data/AGRF_Genotypes_2017_wide.csv")

# load markers information table
Arab_markers <- read_csv("./data/Yasir_marker_export.csv") # %>% .[match(locNames(yasir_genind), .$Locus),]
marker_names <- Arab_markers$Locus
              
SSR_data_2017_pop <- SSR_data_2017_wide  %>% mutate_at(vars(one_of(marker_names)), as.numeric) #%>% 
  # mutate(Genotyped_markers=rowSums(.[marker_names]>0)) # %>% 
  # Change failed and NAs to 0s
# mutate_at(vars(one_of(marker_names)), ~case_when(is.na(.) ~ "0", .=="Failed" ~ "0", TRUE ~ as.character(.))) %>% 
# Add to the database
xlsx::write.xlsx(as.data.frame(SSR_data_2017_pop), "../../data/4_year_basic_complete_file_IB_S.xlsx", 
                 sheetName = "A_rab_SSR_2017_Feb2019", append = TRUE, 
                 row.names = FALSE)
# SSR_data_2017_pop %>% filter(Genotyped_markers==7) # %>% filter(n>1)

#### match haplotype ID with allele combination ####
# load haplotype information and save to file
haplo_table <- read_excel("../../data/Yasir_isolate_haplotype_info.xlsx",  
                          trim_ws = TRUE) %>% mutate(Isolate=sub("^TRTR", "TR", Isolate)) %>% 
  select(Haplotype, Isolate)
Arab_all <- readxl::read_excel("../../data/4_year_basic_complete_file_IB_S.xlsx", col_names = colnames(SSR_data_2017_wide), skip = 1) %>% bind_rows(SSR_data_2017_pop) %>%  mutate(Region=paste0("Reg", Region)) %>% 
  left_join(haplo_table, by=c("Sample_Name"="Isolate")) %>% filter(Year>2010) %>% 
  mutate_if(is.character, ~case_when(grepl("unknown", ., ignore.case = TRUE) ~ NA_character_, TRUE ~ .)) %>% 
  mutate_at(vars(starts_with("Host")), ~sub("^Genesis09[0-9]+$", "Genesis090", 
                                            sub("Genesis ", "Genesis", sub("PBA ", "", .))))# %>% 
  
summary_table <- Arab_all %>% select(Region,Isolate=Sample_Name,Location,State,Host,Year) %>% group_by(Isolate) %>% 
  slice(1) %>% ungroup() %>% group_by(Year) %>% count(Region) %>% 
   spread(Region, n) %>% mutate_all(~replace_na(., 0)) %>%  
  write_xlsx(., "./output/A_rabiei_sampling_collection_2013-2017.xlsx", "collection_sum")



# mutate(Host=sub("PBA ", "", Host)) 
# Check inconsistencies between isolate data (State-based)
# read_excel("../../data/Yasir_isolate_haplotype_info.xlsx",  
#                          trim_ws = TRUE) %>% mutate(Isolate=sub("^TRTR", "TR", Isolate)) %>% 
#   right_join(Arab_pop %>% select(-one_of(marker_names)), by=c("Isolate"="Sample_Name", "Year", "State", "Host", "Region")) %>% 
#   filter(State.x!=State.y | Year.x!=Year.y | Host.x!=Host.y) %>% 
#   select(-one_of(marker_names)) %>% select(order(colnames(.))) %>% 
#   write_csv("./output/A_rabiei_SSR_2017_inconsistent_isolate_info.csv")
# 
# Look for duplicate haplotypes
Arab_all %>% select(Haplotype, one_of(marker_names)) %>% distinct() %>% arrange_at(marker_names) %>% print(n=Inf)
# duplicate samples
Arab_all %>% count(Sample_Name) %>% filter(n>1)
# # Fix inconsistencies manually
# 
# Arab_pop %>% filter(is.na(Host) | is.na(Year))
# load population table
Arab_pop <- Arab_all %>%  mutate(Genotyped_markers=rowSums(.[marker_names]>0), Region=factor(Region)) %>%
  # Filter isolates with missing host and year information and more than 2 missing markers
  dplyr::filter(!is.na(Host), Genotyped_markers>=5, !is.na(Year)) %>%  # , Year==2016
  rename(Ind=Sample_Name) %>% #mutate(Region=paste0("Reg",Region ))  %>% 
   arrange(Ind, desc(Genotyped_markers))# %>% column_to_rownames("Ind")# skip = 2, 

# Add haplotype info and save to file
# haplo_table  %>% filter(!Isolate %in% Arab_pop$Ind) %>% 
#   select(Haplotype) %>% distinct() %>% arrange(Haplotype) %>% print(n=Inf)
Arab_haplos <- haplo_table %>%  inner_join(Arab_pop %>% select(-Haplotype), by=c("Isolate"="Ind")) %>% 
  group_by(Haplotype) %>% # , "Region", "Year"
  slice(1) %>% arrange(desc(Genotyped_markers)) %>% ungroup()
# save to file
xlsx::write.xlsx(as.data.frame(Arab_haplos[c("Haplotype", marker_names)]), 
                 "./output/A_rabiei_SSR_haplotypes.xlsx", sheetName = "common_haplotypes", 
                 row.names = FALSE)
# If there's only one possible allele (other than 0 for missing), replace the missing ones with that allele

set_locus <- function(vec){
  # vec[vec==285] <- 283
  if (length(unique(vec[vec!=0]))==1) return(unique(vec[vec!=0]))
  return(vec)
}


# read_excel("../../data/Yasir_isolate_haplotype_info.xlsx", range = cell_cols(1:4)) %>% 
#   arrange(desc(Haplotype)) %>% slice(1) %>% select(Haplotype)

all_haplos <- Arab_pop[c(marker_names)] %>%  # lapply(., set_locus) %>% as.data.frame() %>% 
  mutate_all(set_locus) %>%
  # mutate(Genotyped_markers=rowSums(.[marker_names]>0)) %>% 
  distinct() %>% full_join(Arab_haplos[c("Haplotype", marker_names)], .) %>% 
  arrange(Haplotype) 
missing_haplos <- all_haplos[is.na(all_haplos$Haplotype),]
last_haplo <- max(as.numeric(str_extract(all_haplos$Haplotype, "\\d+")), na.rm = TRUE)
all_haplos$Haplotype[is.na(all_haplos$Haplotype)] <- paste0("ARH", last_haplo+seq_along(missing_haplos$Haplotype))
  
all_haplos %>% print(n=Inf)
# find duplicates
all_haplos <- all_haplos %>% unite(allele_combined,marker_names) %>% distinct(allele_combined, .keep_all=TRUE) %>% 
  separate(allele_combined, marker_names) %>% mutate_at(vars(one_of(marker_names)), as.numeric)
  # arrange(desc(Genotyped_markers)) %>% 
all_haplos %>%  as.data.frame() %>% 
  xlsx::write.xlsx(., "./output/data/A_rabiei_SSR_haplotypes.xlsx", 
                   sheetName = "all_haplotypes",  row.names = FALSE, append = TRUE)

#### Combine current and previous data with haplotype information #####
Arab_pop <- Arab_pop %>% select(-Haplotype) %>% # unite(allele_combined,marker_names) %>%
  left_join(all_haplos) # %>% unite(allele_combined,marker_names))

# left_join(all_haplos) %>% print(n=Inf)

# Check that there are no missing haplotypes
unique(Arab_pop$Haplotype)


##### QC of data #####
# find duplicate IDs
dup_samples <- Arab_pop %>% count(Ind) %>% filter(n>1)
Arab_pop %>% filter(Ind %in% dup_samples$Ind)  %>% as.data.frame() %>% 
  xlsx::write.xlsx(., "./output/A_rabiei_SSR_2017_duplicate_samples.xlsx", 
                   sheetName = "duplicate_samples", 
                    row.names = FALSE)
  # write_csv("./output/A_rabiei_SSR_2017_duplicate_samples.csv")
# check consistency of Hosts 
unique(Arab_pop$Host)
Arab_pop %>% count(Host) %>% as.data.frame() %>% 
  xlsx::write.xlsx(., "./output/A_rabiei_SSR_2017_duplicate_samples.xlsx", 
                   sheetName = "hosts", 
                   append = TRUE, row.names = FALSE)

# check consistency of Years
Arab_pop %>% count(Year)

# write back to the database
Arab_pop <- Arab_pop %>% # filter() %>% 
  group_by(Ind) %>% arrange(Ind, desc(Genotyped_markers)) %>% slice(1) %>% 
  ungroup() %>% select(-Genotyped_markers) 
Arab_pop %>% as.data.frame() %>% 
  xlsx::write.xlsx(., "./data/5_year_complete_SSR_db.xlsx", sheetName = "Feb_2019", #append = TRUE, 
                   row.names = FALSE)

# check Years
Arab_pop %>% count(Year)
# check regions
Arab_pop %>% count(Region)
# check haplotypes
Arab_pop %>% count(Haplotype)
# chek haplotypes in 2017
Arab_pop %>% filter(Year==2017) %>% count(Haplotype)
all_haplos %>% mutate(genotyped_loci=rowSums(.[marker_names]>0)) %>% arrange(desc(genotyped_loci)) %>% 
  count(genotyped_loci)
  # write_csv("./output/A_rabiei_SSR_2017_hosts.csv")

create_report(Arab_pop %>% mutate(Year=as.character(Year), Haplotype=fct_lump(Haplotype %>% as.factor, 10), Host=fct_lump(Host %>% as.factor, 10)) , output_file = "A_rabiei_SSR_descriptive_report_other_haplo.html", 
              output_dir = "./output",
                config = list("introduce" = list(),"plot_histogram" = list(),
                              "plot_bar" = list()))

#### haplotype distribution in years ####
# plot haplotype distribution per year as a function of region and host
# focus on common haplotypes
common_haplos <- Arab_pop  %>%  count(Haplotype) %>% arrange(desc(n)) %>% slice(1:12)
Arab_pop_common <- Arab_pop %>%
  mutate_at(c("Haplotype", "Host"), ~fct_lump(.x %>% as.factor, n = 10)) # %>%

haplo_freq_by_reg <- Arab_pop_common %>% filter(Year>2010) %>%  group_by(Year, Region, Haplotype) %>% summarise(n=n()) %>% mutate(sum=sum(n), freq=n/sum)

# Plot frequencies by year and region
ggplot(haplo_freq_by_reg, aes(x=Region, y=freq, fill=Haplotype)) + geom_bar(stat = 'identity', width = 0.7) + 
  scale_fill_brewer(palette = "Paired") + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) +
  facet_wrap(~Year) + geom_text(aes(y=0.9, label=sum), size=4.5) + labs(y="Frequency")+
  coord_flip()
ggsave("./plots/A_rabiei_haplotypes_by_region_year_hori.pdf", width = 12, height=7)  

Arab_pop_common %>% # filter(Year>2010)  %>% 
  count(Haplotype) %>% mutate(sum=sum(n), freq=n/sum) %>% arrange(desc(n))

# Plot haplotype frequencies by year 
haplo_freq_by_year <- Arab_pop_common %>% filter(Year>2010) %>%  group_by(Year) %>% count(Haplotype) %>% 
  mutate(sum=sum(n), freq=n/sum)
haplo_labs <- haplo_freq_by_year %>% filter(Haplotype=="ARH01") %>% mutate(label_pos=1-freq, label=sprintf("%.2f%%", freq*100))

ggplot(haplo_freq_by_year, aes(x=Year, y=freq, fill=Haplotype)) + geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_brewer(palette = "Paired") + plot_theme(baseSize = 21) + 
  # theme(strip.text=element_text(size=rel(0.7)),
  #       axis.text.x = element_text(size=rel(0.75))) +
  geom_text(data =haplo_labs, aes(x=Year, y=0.95, label=glue::glue("n={sum}")), size=4.5) +
  geom_text(data =haplo_labs, aes(x=Year, y=label_pos, label=label), size=4.5, nudge_y = 0.065) +
  theme(aspect.ratio = 1/2) + labs(y="Frequency") +
  coord_flip()
# save plot
ggsave("./plots/A_rabiei_haplotypes_year_hori.pdf", width = 10, height=8)  


#### Pathogenicity Analysis ####
# path_levels <- set_names(as.character(c(0,0,1,3)), c("N/A", "Low", "Moderate", "High"))
# diff_set <- c("Seamer", "ICC3996", "Genesis090", "HatTrick")
scoring_table <- readxl::read_excel("../../../Pathology/scoring_table.xlsx", "Path_score") %>% 
  filter(Host!="Seamer") %>% gather(key = "path_levels", value = "score", -Host) 
scoring_set <- unique(scoring_table$Host)
assign_score <- function(path_level, host, scoring_mat=scoring_table){
  scoring_mat %>% filter(Host==host, path_levels %in% path_level) %>% select(score) %>% as.numeric()
}

patho_files <- list.files("../../../Pathology", "A_rabiei_pathogenicity_results.+", 
                          full.names = TRUE) %>% 
  .[!grepl("/~$", ., fixed = TRUE)]
patho_table <- patho_files %>% map_dfr(~readxl::read_excel(., sheet = "patho_table") %>% 
                       mutate_all(as.character) %>% mutate_all(~str_replace_na(., "N/A")) %>% 
                       mutate(Year=as.numeric(str_extract(.x, "20[0-9]{2}"))) %>% 
            set_names(., sub("(.+?) .+", "\\1", colnames(.))))  %>% 
  mutate(Isolate=gsub(" ", "", Isolate),
         Host=sub("^Genesis09[0-9]+$", "Genesis090", 
                  sub("Genesis ", "Genesis", sub("PBA ", "", Host))))

# patho_table$Path_rating <- patho_table %>% select(scoring_set) %>% map(scoring_set, ~)
  
patho_table$Path_rating <- rowSums(map_dfc(scoring_set, 
                                     function(n) map_dbl(patho_table[[n]], ~assign_score(., n))))
patho_table <- patho_table %>% mutate(Path_rating=case_when(#Path_rating>4 ~ 4, 
                                                ICC3996=="High" ~ 5, # double check this with Rebecca
                                      ICC3996!="High" & Genesis090!="High" & HatTrick!="High" ~ 0,
                                                TRUE ~ Path_rating)) %>% 
  mutate(Rating=if_else(Path_rating>0, "High", "Low"), Pathotype=paste0("Group", Path_rating )) %>%    mutate_all(~str_replace_na(., "N/A")) %>% 
  # as.data.frame() %>%
  write_xlsx(., "./output/A_rabiei_pathotypes_2013-2018.xlsx", sheet = "pathotyping")

# patho_table %>% filter(!grepl("N/A|Low", Pathotype), Pathotype!=Path_rating)
# patho_table %>% arrange(desc(Path_rating))
# Pathotype=factor(paste0("Group", Pathotype)), 
#Pathogenicity=if_else(Pathotype==0, "Low", "High"), 
                                         
# Check how many missing from 2017 table
# patho_table %>% filter(Year==2017, Isolate %in% Arab_pop$Ind, Rating=="High")
patho_table %>% count(Year)
# path_levels <- c("Low Pathogenic", paste0("Group", 1:4))
# 
# path_cols <- setNames(brewer_pal(type = "div"), path_levels)
# brewer.pal(5, "Set1")[c(2,3,5,1,4)]
# summarise 
patho_sum <- patho_table %>%  group_by(Year, Pathotype) %>% 
  summarise(n = n()) %>%
  mutate(pop=sum(n), freq = n / pop) %>% 
  arrange(Year, desc(Pathotype)) %>% 
  mutate(label=sprintf("%.1f%%", freq*100), cumsum=cumsum(freq),
         label_pos=freq/2+dplyr::lag(cumsum, default = 0)) %>% as.data.frame() %>% 
  write_xlsx(., excel_file = "./output/A_rabiei_pathotypes_2013-2018.xlsx", sheet = "pathotype_freq", append=TRUE)

patho_labs <- patho_sum %>% filter(Pathotype=="Group0") %>% mutate(label_pos=1-freq)
  # mutate(Pathotype=factor(Pathotype, levels = c("Low Pathogenic", paste0("Group", 1:4))))
 #group_by(Year) %>% 
  # filter(Pathotype %in% paste0("Group", 3:4)) %>% 
  

# statistically check differences in each group
# summarise per year and pathotype
patho_sum %>% filter(Year==2015 | Year==2017) %>% group_by(Pathotype) %>% summarise(prop_test=prop.test(x=n, n=pop, alternative = "two.sided", correct = FALSE)$p.value)



# Plot distribution of pathotype groups
p <- ggplot(patho_sum, aes(x=Year, y=freq, fill=Pathotype))
p + geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_brewer(palette = "Spectral", direction = -1) + plot_theme() +
  scale_x_continuous(breaks= as.integer(unique(patho_sum$Year))) +
  geom_text(data =patho_labs, aes(x=Year, y=0.95, label=glue::glue("n={pop}")), size=4.5) + 
  geom_text(data =patho_labs, aes(x=Year, y=label_pos, label=label), size=4.5, nudge_y = 0.065) +
  theme(aspect.ratio = 1/2) + labs(y="Frequency") +
  coord_flip()

ggsave("./plots/A_rabiei_pathotype_groups_2013-2018_horiz.pdf", width =10, height = 8)

patho_table %>% filter(ICC3996!="Low")
  
diff_hosts <- c(scoring_set, "Seamer")
Patho_host_sum <- patho_table %>% gather(key="Diff_host", value="Diff_rating", diff_hosts) %>% 
  filter(Diff_rating!="N/A") %>%
  group_by(Year, Diff_host) %>% mutate(Diff_rating=str_replace(Diff_rating, "Moderate", "High")) %>%
  count(Diff_rating) %>%  
  mutate(pop=sum(n), freq = n / pop) # %>% ungroup()
  # arrange(Year, desc(Pathotype)) %>% 
  # mutate(label=sprintf("%.1f%%", freq*100), cumsum=cumsum(freq),
  #        label_pos=freq/2+dplyr::lag(cumsum, default = 0)) %>% as.data.frame() %>% 
  # write_xlsx(., excel_file = "./output/A_rabiei_pathotypes_2013-2017.xlsx", sheet = "pathotype_freq", append=TRUE)

ggplot(data=Patho_host_sum %>% filter(Diff_rating=="High"), 
       aes(x=Year, y=freq, group = Diff_host, colour=factor(Diff_host, levels = c("Genesis090", "ICC3996", "HatTrick", "Seamer")))) + geom_line(size=1.5) + labs(y="Frequency of isoaltes capable of causing\n'High' and 'Moderate' disease symptoms", colour="Host") +
  plot_theme(baseSize = 20) + theme(axis.title.y = element_text(size = rel(0.7), face = "plain")) +
  geom_text(data=Patho_host_sum %>% filter(Diff_rating=="High", Year>2015, Year<2018),
            aes(label=sprintf("%.1f%%", freq*100)), size=4, nudge_y = 0.015, nudge_x = -0.15, colour="black") +
  geom_text(data=Patho_host_sum %>% filter(Diff_rating=="High", Year==2018),
            aes(label=sprintf("%.1f%%", freq*100)), size=4, nudge_y = 0, nudge_x = 0.25, colour="black") +
  scale_color_brewer(palette = "Set1")
# save plot
ggsave("./plots/A_rabiei_host_interaction_2013-2018a.pdf", width = 10, height = 7)


#### Haplotype-pathotype analysis ####
path_haplo <- patho_table %>% 
  select(Isolate, Rating, Pathotype, Year) %>% 
  left_join(Arab_pop, by=c("Isolate"="Ind", "Year")) %>% # as.data.frame() %>% 
  write_xlsx(., "./output/A_rabiei_pathotypes_2016-2017.xlsx", "patho_haplo", append=TRUE)
# save inconsistencies to file
patho_table %>% 
  select(Isolate, Rating, Pathotype, Host, Year) %>% 
  inner_join(Arab_pop, by=c("Isolate"="Ind")) %>% 
  filter(Year.x != Year.y | Host.x != Host.y) %>% # as.data.frame() %>% 
  # select(order(colnames(.))) %>% # , "Host", "Year"
  write_xlsx(., "./output/A_rabiei_pathotypes_2016-2017.xlsx", "patho_haplo_inconsist",
                   row.names = FALSE, append = TRUE)

# check if missing isolates were not genotyped completely against the full table
patho_table %>% 
  select(Isolate, Rating, Pathotype, Year) %>% filter(Year==2017) %>% 
  inner_join(SSR_data_2017_pop, by=c("Isolate"="Sample_Name", "Year")) %>% 
  filter(rowSums(.[marker_names]>0)>=5) %>% count(Isolate) %>% filter(n>1)

# inner join of pathology and haplotype data
path_haplo_inner <- patho_table %>% 
  select(Isolate, Rating, Pathotype, Year) %>% 
  inner_join(Arab_pop, by=c("Isolate"="Ind", "Year"))

common_haplotypes <- path_haplo_inner %>% count(Haplotype) %>% filter(n>1)
path_haplo_inner %>% count(Haplotype) %>% arrange(desc(n)) %>% print(n=Inf)
common_hosts <- path_haplo_inner %>% count(Host) %>% filter(n>1)
# combine rare hosts and haplotypes 
path_haplo_combined <- path_haplo_inner  %>% select(-marker_names) %>% 
  mutate(Haplotype=fct_lump(Haplotype %>% as.factor, n = 11)) %>% 
  mutate(Host=fct_lump(Host %>% as.factor, n = nrow(common_haplotypes)), Rating=factor(Rating, levels=c("Low", "High")))
haplo_cols <- setNames(RColorBrewer::brewer.pal(12, "Paired"  ), levels(path_haplo_combined$Haplotype))

# summarise and plot per haplotype Pathotype in group
path_haplo %>% count(Year)
patho_table %>% count(Year)
patho_haplo_sum <- path_haplo_combined %>% group_by(Year, Pathotype) %>% count(Haplotype) %>% 
   mutate(sum=sum(n), freq=n/sum)

ggplot(patho_haplo_sum, aes(x=Pathotype, y=freq, fill=Haplotype)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = haplo_cols[patho_haplo_sum$Haplotype]) + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) + coord_flip() +
  facet_wrap(~Year) + labs(y="Frequency") + 
  geom_text(data=patho_haplo_sum %>% filter(Haplotype=="ARH01"), aes(y=0.9, label=sum), size=4.5) 
  
# save plot
ggsave("./plots/A_rabiei_haplotypes_by_year_groups_hori.pdf", width = 12, height=7)

# summarise and plot per haplotype Rating in group 
patho_haplo_rating <- path_haplo_combined %>% group_by(Year, Rating) %>% count(Haplotype) %>% 
  mutate(sum=sum(n), freq=n/sum)
# summarise per year and pathotype
patho_haplo_rating %>% filter(Haplotype=="ARH01") %>% group_by(Year, Haplotype) %>% summarise(prop_test=prop.test(x=n, n=sum, alternative = "less", correct = FALSE)$p.value)

ggplot(patho_haplo_rating, aes(x=Rating, y=freq, fill=Haplotype)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = haplo_cols[patho_haplo_sum$Haplotype]) + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) + coord_flip() +
  facet_wrap(~Year) + labs(y="Frequency", x="Pathogenicity Rating") + theme(aspect.ratio = 1/1.25) +
  geom_text(data=patho_haplo_rating %>% filter(Haplotype=="ARH01"), aes(y=0.9, label=sum), size=4.5) +
  geom_text(data=patho_haplo_rating %>% filter(Haplotype=="ARH01"), aes(y=1-freq, label=percent(freq)), 
            size=4.5, nudge_y = 0.12)

# save plot
ggsave("./plots/A_rabiei_haplotypes_by_year_rating_hori.pdf", width = 12, height=7)  


# summarise and plot pathotype haplotypes per host
Hosts <- c("Genesis090", "Seamer", "HatTrick")
host_haplo_sum <- path_haplo_combined %>% group_by(Host, Pathotype) %>% count(Haplotype) %>% 
  filter(Host %in% Hosts) %>% mutate(sum=sum(n), freq=n/sum) %>% ungroup() %>% mutate(Host=factor(Host, levels = Hosts)) 

ggplot(host_haplo_sum, aes(x=Pathotype, y=freq, fill=Haplotype)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = haplo_cols[host_haplo_sum$Haplotype]) + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) + coord_flip() + #theme(aspect.ratio = 1/1.25)+
  facet_wrap(~Host, ncol = 2, drop = TRUE) + labs(y="Frequency") + 
  geom_text(data=host_haplo_sum %>% filter(Haplotype=="ARH01"), aes(y=0.9, label=sum), size=4.5) 

# save plot
ggsave("./plots/A_rabiei_haplotypes_by_host_groups_hori.pdf", width = 12, height=7)

# summarise and plot rating haplotypes per host
Hosts <- c("Genesis090", "Seamer", "HatTrick")
host_haplo_rating <- path_haplo_combined %>% group_by(Host, Rating) %>% count(Haplotype) %>% 
  filter(Host %in% Hosts) %>% mutate(sum=sum(n), freq=n/sum) %>% ungroup() %>% mutate(Host=factor(Host, levels = Hosts)) 
# check statistical difference between high and low in each host
host_haplo_rating %>% filter(Haplotype=="ARH01") %>% group_by(Host, Haplotype) %>% summarise(prop_test=prop.test(x=n, n=sum, alternative = "less", correct = FALSE)$p.value)
# check statistics Genesis090 and Seamer
host_haplo_rating %>% filter(Haplotype=="ARH01", grepl("Gene|Seamer", Host)) %>% group_by(Rating, Haplotype) %>% summarise(prop_test=prop.test(x=n, n=sum, alternative = "greater", correct = FALSE)$p.value)
# check statistics HatTrick and Seamer
host_haplo_rating %>% filter(Haplotype=="ARH01", grepl("Hat|Seamer", Host)) %>% group_by(Rating, Haplotype) %>% summarise(prop_test=prop.test(x=n, n=sum, alternative = "greater", correct = FALSE)$p.value)

# plot rating haplotypes per host
ggplot(host_haplo_rating, aes(x=Rating, y=freq, fill=Haplotype)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = haplo_cols[host_haplo_rating$Haplotype]) + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) + coord_flip() + #theme(aspect.ratio = 1/1.25)+
  facet_wrap(~Host, ncol = 2, drop = TRUE) + labs(y="Frequency", x="Pathogenicity Rating") + 
  geom_text(data=host_haplo_rating %>% filter(Haplotype=="ARH01"), aes(y=0.9, label=sum), size=4.5) +
  geom_text(data=host_haplo_rating %>% filter(Haplotype=="ARH01"), aes(y=1-freq, label=percent(freq)), 
            size=4.5, nudge_y = 0.08)

# save plot
ggsave("./plots/A_rabiei_haplotypes_by_host_rating_hori.pdf", width = 12, height=7)  

# summarise and plot per host (pie chart)
Hosts <- c("Genesis090", "Seamer", "HatTrick")
host_haplo_sum <- path_haplo_combined %>% group_by(Host) %>% count(Pathotype) %>% 
  filter(Host %in% Hosts) %>% mutate(sum=sum(n), freq=n/sum) %>% ungroup() %>% mutate(Host=factor(Host, levels = Hosts)) 

ggplot(host_haplo_sum, aes(x=Pathotype, y=freq, fill=Haplotype)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  scale_fill_manual(values = haplo_cols[host_haplo_sum$Haplotype]) + plot_theme(baseSize = 21) + 
  theme(strip.text=element_text(size=rel(0.7)),
        axis.text.x = element_text(size=rel(0.75))) + coord_flip() + #theme(aspect.ratio = 1/1.25)+
  facet_wrap(~Host, ncol = 2, drop = TRUE) + labs(y="Frequency") + 
  geom_text(data=host_haplo_sum %>% filter(Haplotype=="ARH01"), aes(y=0.9, label=sum), size=4.5) 

# save plot
ggsave("./plots/A_rabiei_haplotypes_by_host_groups_hori.pdf", width = 12, height=7)  

##### Donut around pie ############
# from http://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot
givemedonutsorgivemedeath(family_summary, "phylum", "family", "log_count",inner_labs = T, legend_offset = c(0.025,-0.09))


#### Summary table ####
patho_table %>% group_by(Year) %>% count(State)
Arab_pop_common %>% group_by(Year) %>% count(Haplotype) %>% mutate(sum=sum(n), freq=n/sum)


#### poppr diversity analysis #####
# colnames(Arab_pop)[1:2] <- c("Region", "Ind")
Arab_genind <- df2genind(Arab_pop[marker_names], ploidy = 1, NA.char = 0, # ncode = 3,
                          ind.names = Arab_pop$Ind,
                          pop =  Arab_pop$Region, 
                         strata = Arab_pop[!colnames(Arab_pop) %in% c("Ind", marker_names)])#), loc.names = Arab_markers$Locus) #, 
                          # strata = )
# Need to clean some individuals with missing information


clean_genind <- Arab_genind %>% missingno("genotype", cutoff=0.51)
clean_Ind <- which(Arab_pop$Ind %in% indNames(clean_genind))
# Add host information
# 
# strata_df <- Arab_pop[clean_Ind,!colnames(Arab_pop) %in% c("Ind", Arab_markers$Locus)]
# strata(clean_genind) <- strata_df

LogMsg(sprintf("Microsatellite data imported into 'genind' object.\nOverall, %d individuals from %d regions and across %d years were imported, with %d loci.",
            nInd(clean_genind), nPop(clean_genind), length(levels(strata(clean_genind)$Year)), nLoc(clean_genind)))

# names(Arab_pops)

#### Data QC and EDA ####
pdf("./plots/Genotype_accumulation_curve.pdf", width = 8, height = 6) 
genotype_curve(clean_genind, sample = 1000, quiet = TRUE)
dev.off()
# Look at locus details
(locus_stats <- locus_table(clean_genind))

alleles_per_region <- Arab_pop %>% gather(key="locus", value = "Allele", marker_names) %>% 
  filter(Allele!=0) %>% group_by(Region, locus) %>% count(Allele) %>% count(locus) %>% spread(Region, n)

locus_stats %>% as.data.frame() %>% rownames_to_column("locus") %>% filter(locus!="mean") %>% 
  left_join(alleles_per_region) %>% 
  mutate(allele_sizes=alleles(clean_genind) %>% map(~paste(sort(.), collapse = ","))) %>% 
  xlsx::write.xlsx(., "./output/A_rabiei_SSR_poppr.xlsx", sheetName = "locus_stats", row.names = FALSE)

Arab_pop %>% filter(Ind=="14HOR017")
# isolates showing "new" alleles
Arab_pop %>% filter(`ME14-1-56`=="381" | ArH05T=="236")

#### Population diversity ####
Arab_pop_2017 <- Arab_pop %>% filter(Year==2017)
pop_gen_2017 <-  df2genind(Arab_pop_2017[marker_names], ploidy = 1, NA.char = 0, # ncode = 3,
                                                              ind.names = Arab_pop_2017$Ind,
                                                              pop = Arab_pop_2017$Region, 
                                                              strata = Arab_pop_2017[c("Region", "State", "Host", "Haplotype")])
poppr(pop_gen_2017) %>% mutate(lambda_corr=N/(N - 1) * lambda) %>% arrange(Pop) %>% 
  write_xlsx(., "./output/A_rabiei_SSR_poppr.xlsx", "by_region_2017", append = TRUE)

# Generate population genetics metrices
region_pop_gen <- poppr(clean_genind) %>% mutate(lambda_corr=N/(N - 1) * lambda) %>% arrange(Pop) %>% 
  write_xlsx(., "./output/A_rabiei_SSR_poppr.xlsx", "by_region_all_years", append = TRUE)

Arab_by_years <- clean_genind
pop(Arab_by_years) <- strata(clean_genind)$Year
Arab_by_years <- Arab_by_years[!is.na(pop(Arab_by_years))]

# Generate population genetics metrices
year_pop_gen <- poppr(Arab_by_years) %>% mutate(lambda_corr=N/(N - 1) * lambda) %>% arrange(Pop) %>% 
  write_xlsx(., "./output/A_rabiei_SSR_poppr.xlsx", "by_year_all_regions", append = TRUE)



# Recombine just the wanted populations
# selected_pops <- repool(Arab_pops[[Pops[1]]], Arab_pops[[Pops[2]]])


# Export to Genalex
genind2genalex(clean_genind, filename = "./output/A_rabiei_SSR_2017_by_Host_genalex.csv", pop=strata(clean_genind)$Host)
genind2genalex(clean_genind, filename = "./output/A_rabiei_SSR_2017_by_Year_genalex.csv", pop=strata(clean_genind)$Year)
genind2genalex(clean_genind, filename = "./output/A_rabiei_SSR_2017_by_Region_genalex.csv", pop=pop(clean_genind))

#### plot haplotype Networks ####
# define populations to use
# sub_Pops <- c("HatTrick", "Genesis090")
# separate by host
cols <- setNames(rainbow(nPop(clean_genind)), popNames(clean_genind))
hosts_cols <- setNames(rainbow(length(Hosts)), Hosts)
# clean_pops <- seppop(clean_genind)
# Recombine just the wanted populations
# selected_pops <- clean_pops[[sub_Pops[1]]]
# Plot msn for each population
for (h in Hosts){
  # h <- Hosts[1]
  pe_data <- clean_genind[clean_genind@strata$Host==h,] # %>%  missingno("genotype", cutoff=0.2)
  # Calculate the MSN
  pe_data.msn <- bruvo.msn(pe_data, replen = Arab_markers$Repeat_size,  showplot = FALSE)
  # Visualize the network
  set.seed(120)
  pdf(filedate(sprintf("scaled_bruvo.msn_regions_plot_%s_only", h), ".pdf", "plots"), width = 10, height = 7)
  plot_poppr_msn(pe_data, pe_data.msn, inds = "none", palette = cols,  nodescale = 10, 
                 main=sprintf("bruvo.msn plot of haplotype distribution in growing regions (%s only)", h))# , nodescale = 5 nodebase = 1.25,
  dev.off()
}


# For all data combined

pe_data <- clean_genind # %>%  missingno("genotype", cutoff=0.2)
# Calculate the MSN
pe_data.msn <- bruvo.msn(pe_data, replen = Arab_markers$Repeat_size,  showplot = FALSE)
# Visualize the network
set.seed(120)
pdf(filedate(sprintf("scaled_bruvo.msn_regions_plot_all_hosts"), ".pdf", "plots"), width = 10, height = 7)
plot_poppr_msn(pe_data, pe_data.msn, inds = "none", palette = cols,  nodescale = 10,
               main="bruvo.msn plot of haplotype distribution in growing regions (all hosts)")#, nodebase = 1.25) #  
dev.off()

# Plot recent year
# Calculate the MSN
pe_data_2017.msn <- bruvo.msn(pop_gen_2017, replen = Arab_markers$Repeat_size,  showplot = FALSE)
# Visualize the network
set.seed(120)
pdf(filedate(sprintf("scaled_bruvo.msn_regions_plot_all_hosts_2017"), ".pdf", "plots"), width = 10, height = 7)
plot_poppr_msn(pop_gen_2017, pe_data_2017.msn, inds = "none", palette = cols,  nodescale = 10,
               main="bruvo.msn plot of haplotype distribution in growing regions (all hosts)")#, nodebase = 1.25) #  
dev.off()


##### Plot dendrogram  ######

# for (p in Pops){
#   p=Pops[1]
  pe_data <- selected_pops %>%
    missingno("genotype", cutoff=0)
  # create a tree using Bruvo's Distance (for microsat markers)
  bruvo.boot(pe_data, replen = Arab_markers$Repeat_size, sample=100, cutoff=80)
  # Import the data as clone data and removing genotypes with 20% missing data
  clone_data <- Arab_pops[[p]] %>%
    missingno("genotype", cutoff=0.2)
 # inds_20 <- sample(nInd(clone_data), 20)
  pop_dist <- aboot(pe_data, dist = provesti.dist, sample = 200,  cutoff = 50) # [inds_20]
  theTree <- pop_dist$edge %>%
    ape::nj() %>%    # calculate neighbor-joining tree
    ape::ladderize() # organize branches by clade
  plot(theTree)
  add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

  # try to paint each population
  # cols <- rainbow()

  ape::plot.phylo(clone_dist, cex = 0.8, font = 1, adj = 0, tip.color = cols[levels(clone_data$pop)],
             label.offset = 0.0125)
  nodelabels(MXtree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
             font = 3, xpd = TRUE)
  axisPhylo(3)

# }


  
