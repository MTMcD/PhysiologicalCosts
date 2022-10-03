# PhysiologicalCosts
Analysis of a field experiment in Barn Swallows with simultaneous brood size manipulation and GPS tagging

## Data files
BSM_FemalePhys_wide.csv: main data file, containing physiological and reproductive information about birds used in the experiment, arranged so that each row corresponds to one individual

BSM_FemalePhys.csv: a subset of physiological and reproductive information about birds used in the experiment, arranged so that each row corresponds to a capture event (so that each individual may have multiple rows)

SiteKey.csv: a list of breeding sites and corresponding grouping regions

filteredGPS19-20.csv: processed GPS data containing data on range size for each tagged individual

## Scripts
BSM_phys_analysis.R: processes field data, and analyzes physiological and reproductive outcomes associated with brood size manipulation and GPS tagging

rawGPS.R: takes a folder of raw GPS data files (such as that stored in MoveBank) and calculates minimum convex polygons for range size at particular chick ages

BSM_GPS_analysis.R: processes field data, and analyzed physiological and reproductive outcomes associated with brood size manipulation and range size

model_output.R: script tidies model output (lmer or glmer) and produces tables and figures with model estimates for publication
