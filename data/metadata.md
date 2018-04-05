# Datasets for fish fecundity and egg sizes  
## data/fecundityFemaleSize.csv  
### Variables descriptions  
*Species:* Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*spawningMode:* Type of spawning as defined in Kasimatis & Riginos 2016 Coral Reefs. One of three types: scatterer, demersal or pelagic;  
*FemaleMass_g:* Female body mass (actual individual or average across a sample of individuals) in grams from which the eggs have been measured;  
*FemaleSize_mm:* Female body size (actual individual or average across a sample of individuals) in millimetres from which the eggs have been measured;  
*sizeType:* Type of female body size, including SL (standard length), FL (fork length), and TL (total length);  
*Fecundity_nOfEggs_per_female:* Total female fecundity (either sampled or estimated) reported in original publication;  
*typeOfWork_Fecundity:* Type of measurement: "Field" means that observations were immediately done from field-collected material. "Field/Lab" means that field-collected material was kept for a determined period (indicated in the respective time column in months) in the lab prior to observations. Blank values are "Field";  
*timeToFecundityMeasurement_months:* For the "Field/Lab" values in the typeOfWork_Fecundity column, indicates the number of months between specimen collection and trait measurement;  
*Location:* Place of collection;  
*Date:* Date of collection (year and month);  
*SampleSize:* Number of eggs measured;  
*Fertilized:* Logical, TRUE if eggs are fertilized, and FALSE if they are not;  
*Live:* Logical, TRUE if eggs are live, and FALSE if they are preserved;  
*Latitude:* Latitude of place of collection in decimals;  
*Longitude:* Longitude of place of collection in decimals;  
*ReferencePdf:* .pdf file from which the data has been collected (contact the authors in case you want to access some or all of these);  
*lnFemaleSize:* natural logarithm of FemaleSize_mm;  
*lnFecundity:* natural logarithm of Fecundity_nOfEggs_per_female;  
*absLatitude:* absolute Latitude.  
  
## data/eggSizeFemaleSize.csv  
### Variables descriptions  
*Species:* Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*spawningMode:* Type of spawning as defined in Kasimatis & Riginos 2016 Coral Reefs. One of three types: scatterer, demersal or pelagic;  
*FemaleMass_g:* Female body mass (actual individual or average across a sample of individuals) in grams from which the eggs have been measured;  
*FemaleSize_mm:* Female body size (actual individual or average across a sample of individuals) in millimetres from which the eggs have been measured;  
*sizeType:* Type of female body size, including SL (standard length), FL (fork length), and TL (total length);  
*eggSize_mm:* Egg diameter in millimetres which generally represent average across multiple measurements. For ellipsoid type eggs, diameter has been obtained by first estimating the volume of an ellipse, and then calculating the diameter of an spherical egg with equivalent volume;  
*typeOfWork_eggSize:* Type of measurement: "Field" means that observations were immediately done from field-collected material. "Field/Lab" means that field-collected material was kept for a determined period (indicated in the respective time column in months) in the lab prior to observations. Blank values are "Field";  
*timeToEggMeasurement_months:* For the "Field/Lab" values in the typeOfWork_eggSize column, indicates the number of months between specimen collection and trait measurement;  
*Location:* Place of collection;  
*Date:* Date of collection (year and month);  
*SampleSize:* Number of eggs measured;  
*Fertilized:* Logical, TRUE if eggs are fertilized, and FALSE if they are not;  
*Live:* Logical, TRUE if eggs are live, and FALSE if they are preserved;  
*Latitude:* Latitude of place of collection in decimals;  
*Longitude:* Longitude of place of collection in decimals;  
*ReferencePdf:* .pdf file from which the data has been collected (contact the authors in case you want to access some or all of these);  
*lnFemaleSize:* natural logarithm of FemaleSize_mm;  
*lnEggSize:* natural logarithm of eggSize_mm;  
*absLatitude:* absolute Latitude.  
  
## data/fecundityEggSizeFemaleSize.csv  
### Variables descriptions  
*Species:* Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*spawningMode:* Type of spawning as defined in Kasimatis & Riginos 2016 Coral Reefs. One of three types: scatterer, demersal or pelagic;  
*FemaleMass_g:* Female body mass (actual individual or average across a sample of individuals) in grams from which the eggs have been measured;  
*FemaleSize_mm:* Female body size (actual individual or average across a sample of individuals) in millimetres from which the eggs have been measured;  
*Fecundity_nOfEggs_per_female:* Total female fecundity (either sampled or estimated) reported in original publication;  
*sizeType:* Type of female body size, including SL (standard length), FL (fork length), and TL (total length);  
*eggSize_mm:* Egg diameter in millimetres which generally represent average across multiple measurements. For ellipsoid type eggs, diameter has been obtained by first estimating the volume of an ellipse, and then calculating the diameter of an spherical egg with equivalent volume;  
*typeOfWork_eggSize:* Type of measurement: "Field" means that observations were immediately done from field-collected material. "Field/Lab" means that field-collected material was kept for a determined period (indicated in the respective time column in months) in the lab prior to observations. Blank values are "Field";  
*timeToEggMeasurement_months:* For the "Field/Lab" values in the typeOfWork_eggSize column, indicates the number of months between specimen collection and trait measurement;  
*Location:* Place of collection;  
*Date:* Date of collection (year and month);  
*SampleSize:* Number of eggs measured;  
*Fertilized:* Logical, TRUE if eggs are fertilized, and FALSE if they are not;  
*Live:* Logical, TRUE if eggs are live, and FALSE if they are preserved;  
*Latitude:* Latitude of place of collection in decimals;  
*Longitude:* Longitude of place of collection in decimals;  
*ReferencePdf:* .pdf file from which the data has been collected (contact the authors in case you want to access some or all of these);  
*lnFemaleSize:* natural logarithm of FemaleSize_mm;  
*lnEggSize:* natural logarithm of eggSize_mm;  
*lnFecundity:* natural logarithm of Fecundity_nOfEggs_per_female;  
*absLatitude:* absolute Latitude.  

# Dataset for fish egg energy content  
This dataset was collected and compiled by D R Roberston over the past 3 decades. Check methods described in Robertson & Collin (2015) Front Ecol Evol, 2 | doi: 10.3389/fevo.2014.00084 . The file `data/eggEnergyComplete.csv` contains data for multiple individuals within each species.  

## data/eggEnergyComplete.csv  
### Variables descriptions  
*date_day_month_year:* sample processing date, in day / month / year format;  
*family:* Family name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*Species:* Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*volume_mm3:* Estimated egg volume in mm^3;  
*eggDryWeight_mg:* Egg dry weight in mg;  
*energy_j_per_mg:* dry weight energy density of eggs (J / mg) in that clutch;  
*country:* Country where data were collected;  
*ocean:* ocean basin to where data were collected;  
*dehydrationMethod:* either freeze dried or oven dried;  
*obs:* additional observation with respect to handling of samples.  
  
# data/eggDryWeight.csv  
### Variables descriptions  
*Species*: Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*eggDiameter_mm*: Egg diameter in millimetres;  
*eggDryWeight_mg*: Egg dry weight in milligrams;  
*method*: drying method;  
*ReferencePdf*: .pdf file from which the data has been collected (contact the authors in case you want to access some or all of these);  
*obs*: any general observation pertaining to the data;  
*dataSource*: Indicates the precise location in the original reference from which the data was extracted;  
*Year*: Date of collection (year and month);  
*Location*: Place of collection;  
*Latitude:* Latitude of place of collection in decimals;  
*Longitude:* Longitude of place of collection in decimals.  
  
# Datasets used in Table S8  
## data/trioData*  
### General description  
Datasets in this folder come from 3 studies which contain all four variables of interest (female size, fecundity, egg-volume, and egg-energy) where all four came from the same study and population. Variable names are self-explanatory, and results of mass-scaling of reproductive output are detailed in Table S8.  
  
# Auxiliary Datasets  
## data/otlSppSubstitutions.csv  
### General description  
This dataset is used to substitute original OTT IDs in the Actinopterygii downloaded from the Open Tree of Life with actual species names in the datasets above for which we use phylogenetic hierarchical analyses.  
  
### Variables descriptions  
*Species:*  Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*ottFrame:*  The string of species name that actually matches an existing OTT ID in the Actinopterygii tree from the Open Tree of Life;  
*obs:*  Indicates whether the species name in column ottFrame is a misnomer relative to the Catalog of Fishes.  
  
## data/spawningMode.csv  
### General description  
This dataset is contains the taxonomic information and spawning mode classification for all species analysed in this paper (contained in datasets above).

### Variables descriptions  
*order:*  Taxonomical Order as extracted using R package rfishbase;  
*family:*  Taxonomical Family as extracted using R package rfishbase;  
*genus:*  Genus name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*species:*  Species name following the taxonomy provided in the online [Catalog of Fishes](http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp);  
*spawningModeDetailed:* Type of spawning mode (pelagic, demersal, mouth brooder, scatterer, internal brooding, pouch brooder);  
*spawningModeCode:*  Acronyms for column spawningModeDetailed;  
*source:*  Reference from which the spawning information was obtained.  
