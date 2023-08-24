# TURRIS: an open source database and Python tools to compute the fundamental frequency of historic masonry towers
 
***Authors: Arnaud Montabert*, Cédric Giry*, Claire Limoge Schraen, Jade Lépine, Clarisse Choueiri, E. Diego Mercerat, and Plilippe Guéguen***
*Corresponding authors: arnaud.montabert@ens-paris-saclay.fr, cedric.giry@ens-paris-saclay.fr*

## Scientific motivation
Cultural Heritage buildings are complex but of inestimable value. They requires the synergy of all involved communities. We assembled a large database of historical masonry towers through an extensive literature review in the continuity of previous valuable works (). The database collects information about the towers features and the measured fundamental frequency. Moreover, Python scripts are developped to evaluate their fundamental frequency based on empirical, physics-based formumlation, and also a Rayleigh-Ritz approach. The database and the code are available. We highly encourage the community to provide us the characteristics of new instrumented masonry towers to update the **TURRIS** database. Please contact the corresponding authors. 

## Description of the TURRIS database
**T**owers feat**UR**es & f**R**equenc**I**es databa**S**e (TURRIS) is an open database collecting information of 244 historical masonry towers and information about the associated dynamic identification. Each tower is described by 20 parameters (identification, building name, country, town, latitude, longitude, height, effective height, shape of the section, regularity of the elevation, breadth, length, the minimum and maximum wall thickness, its relation with adjacent building, age, bells mass, Young modulus, mass density, Poisson ratio). Each Operational Modal Analysis survey is characterized by parameters (fundamental frequency, minimum, maximum, and standard deviation of the fundamental frequencies for long-term analysis, the nature of the mode shape, duration of the records, sampling rate, and identification technique). The 244 historical masonry towers are located worldwide (mainly in Europe). The TURRIS database allows studying the relationship between geometrical, material parameters, and fundamental frequency. The details of each field of the TURRIS database are given in the following. 

- Id : the identifier for each study 
- database : the collection where the study was taken into account 
- references : the reference of the scientific study. The reference id maybe used with the associated bibtex file., 
- building_name : the name of the historical building 
- country, town, latitude, longitude : the location of the building 
- input : the origin of the sollicitation 
- f0 [Hz] : the measured fundamental frequency 
- H [m] : the height of the tower 
- Heff  [m]: the effective height of the tower computed as the difference between the absolute height of the building and the height of the any adjacent building
- shape : the shape of the section (SQ : square, REC : rectangular, CIR : circular) 
- elevation_regularity : the regularity of section from bottom to the top 
- width [m]: the shorter dimension of the section 
- length [m]: the larger dimension of the section 
- min_wall_thickness [m]: the minimum thickness of the wall of the tower 
- max_wall_thickness [m]: the maximum thickness of the wall of the tower 
- relation : the type of connection between the tower and adjacent building (Isolated or bounded) 
- period : the age of construction 
- bells [kg]: the mass of the bell system 
- E [GPa]: the Young modulus 
- density [kg.m-3]: the mass density 
- Poisson_ratio : the Poisson ratio 
- duration [minute]: the duration of the records 
- sampling [Hz] : the sampling rate of the records 
- info : any information about potential damage or retrofitting actions 
- OMA_technique : the technique used to identify modal parameters 
- min_f0 [Hz]: the minimum measured fundamental frequency measured over time  
- max_f0  [Hz]: the maximum measured fundamental frequency measured over time 
- std_f0 [[Hz]: the standard deviation of fundamental frequency measured over time 

 
