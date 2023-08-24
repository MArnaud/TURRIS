# TURRIS: an open source database and Python tools to compute the fundamental frequency of historic masonry towers
 
***Authors: Arnaud Montabert*, Cédric Giry*, Claire Limoge Schraen, Jade Lépine, Clarisse Choueiri, E. Diego Mercerat, and Plilippe Guéguen***

*Corresponding authors: arnaud.montabert@ens-paris-saclay.fr, cedric.giry@ens-paris-saclay.fr*

***How to cite this work*** Montabert, A.; Giry, C.; Limoge Schraen, C.; Lépine, J.; Choueiri, C. Mercerat, E. D.; Guéguen, P. An Open Database to Evaluate the Fundamental Frequency of Historical Masonry Towers through Empirical and Physics-based Formulations. Buildings 2023

![tower_model](./figure/timo_SAM.png)

## Scientific motivation
Cultural Heritage buildings are complex but of inestimable value. They requires the synergy of all involved communities. We assembled a large database of historical masonry towers through an extensive literature review in the continuity of previous valuable works ([1-11]). The database collects information about the towers features and the measured fundamental frequency. Moreover, Python scripts are developped to evaluate their fundamental frequency based on empirical, physics-based formumlation, and also a Rayleigh-Ritz approach. The database and the code are available. We highly encourage the community to provide us the characteristics of new instrumented masonry towers to update the **TURRIS** database. Please contact the corresponding authors. 

## Description of the TURRIS database
**T**owers feat**UR**es & f**R**equenc**I**es databa**S**e (TURRIS) is an open database collecting information of 244 historical masonry towers and information about the associated dynamic identification. Each tower is described by 20 parameters (identification, building name, country, town, latitude, longitude, height, effective height, shape of the section, regularity of the elevation, breadth, length, the minimum and maximum wall thickness, its relation with adjacent building, age, bells mass, Young modulus, mass density, Poisson ratio). Each Operational Modal Analysis survey is characterized by parameters (fundamental frequency, minimum, maximum, and standard deviation of the fundamental frequencies for long-term analysis, the nature of the mode shape, duration of the records, sampling rate, and identification technique). The 244 historical masonry towers are located worldwide (mainly in Europe). The TURRIS database allows studying the relationship between geometrical, material parameters, and fundamental frequency. The details of each field of the TURRIS database are given in the following. 

- *Id* : the identifier for each study 
- *database* : the collection where the study was taken into account 
- *references* : the reference of the scientific study. The reference id maybe used with the associated bibtex file., 
- *building_name* : the name of the historical building 
- *country*, *town*, *latitude*, *longitude* : the location of the building 
- *input* : the origin of the sollicitation 
- *f0* [Hz] : the measured fundamental frequency 
- *H* [m] : the height of the tower 
- *Heff*  [m]: the effective height of the tower computed as the difference between the absolute height of the building and the height of the any adjacent building
- *shape* : the shape of the section (SQ : square, REC : rectangular, CIR : circular) 
- *elevation_regularity* : the regularity of section from bottom to the top 
- *width* [m]: the shorter dimension of the section 
- *length* [m]: the larger dimension of the section 
- *min_wall_thickness* [m]: the minimum thickness of the wall of the tower 
- *max_wall_thickness* [m]: the maximum thickness of the wall of the tower 
- *relation* : the type of connection between the tower and adjacent building (Isolated or bounded) 
- *period* : the age of construction 
- *bells* [kg]: the mass of the bell system 
- *E* [GPa]: the Young modulus 
- *density* [kg.m-3]: the mass density 
- *Poisson_ratio* : the Poisson ratio 
- *duration* [minute]: the duration of the records 
- *sampling* [Hz] : the sampling rate of the records 
- *info* : any information about potential damage or retrofitting actions 
- *OMA_technique* : the technique used to identify modal parameters 
- *min_f0* [Hz]: the minimum measured fundamental frequency measured over time  
- *max_f0*  [Hz]: the maximum measured fundamental frequency measured over time 
- *std_f0* [[Hz]: the standard deviation of fundamental frequency measured over time

## Code description
### Empirical and physics-based formulation to evaluate the fundamental frequency of a tower
The codes to compute the fundamental frequency using empirical and/or physics based formuations are provided in the script *EmpiricalPhysicsBasedModel.py*. The jupyter notebook entitled Empirical_Physical_Relation.ipynb as been written to provide theoritical background and an example of use.

### Evaluating the fundamental frequency using a Rayleigh-Ritz approach
A Rayleigh-Ritz approach based on a Timoshenko formulation is proposed in the python script *Beam.py*. The notebook *Rayleigh_Ritz.ipynb* describe theoritical background, and a description of the code.

## Acknowledgements
The authors kindly acknowledge the institutions and researchers who provided
additional and helpful information of masonry towers features: the municipality of Montboucher-sur-Jabron, Fernando Lopez-Caballero for the fruitful discussion and the advice regarding the sensitivity analysis, Clotilde Chambreuil and Héloïse Rostagni for their help in completing the database. 
This authors wish to express their most grateful thanks to the French National Research Agency (ANR) for the funding of the ACROSS project (ANR-20-CE03–0003) by which a part of this study has been carried out.

## References
**[1]** Lund, J.; Selby, A.; Wilson, J. The dynamics of bell towers-a survey in northeast England. WIT Transactions on The Built Environment 1995, 17, 8.

**[2]** Schmidt, T. Dynamic behaviour of twin bell towers. In Proceedings of the Proceedings of the 2nd international operational modal analysis conference, 2007.

**[3]** Schmidt, T. Dynamic behaviour of twin bell towers. In Proceedings of the Proceedings of the 2nd international operational modal analysis conference, 2007.

**[4]** Rainieri, C.; Fabbrocino, G. Output-only modal identification for prediction of the elastic period of masonry towers. In Proceedings of the Proc. of the 4th International Operational Modal Analysis Conference, 2011.

**[5]** Limoge, C. Méthode de diagnostic à grande échelle de la vulnérabilité sismique des Monuments Historiques: Chapelles et églises baroques des hautes vallées de Savoie: Large-scale seismic vulnerability assesment method for the masonry architectural heritage: baroque chapels and churches of the French Savoye. PhD thesis, Université Paris-Saclay (ComUE), 2016.

**[6]** Ziegler, A. Dynamik der Glockentürme. In Bauwerksdynamik und Erschütterungsmessungen; Springer Fachmedien Wiesbaden, 2017; pp. 153–165. https://doi.org/10.1007/978-3-658-16054-8_6.

**[7]** Ruiz-Jaramillo, J.; Montiel-Vega, L.; García-Pulido, L.J.; Muñoz-González, C.; Blanca-Hoyos, Á. Ambient Vibration as a Basis for Determining the Structural Behaviour of Watchtowers against Horizontal Loads in Southeast Spain. Applied Sciences 2020, 10, 6114. https://doi.org/10.3390/app10176114.

**[8]** Mercerat, D.; Montabert, A.; Giry, C.; Lancieri, M.; Arrighetti, A. Operational Modal Analysis of five historical bell towers in the Mugello basin (Tuscany, Italy). In Proceedings of the Proceedings of the 3 rd European conference on earthquake engineering & seismology (3ECEES), 2022, pp. 4106–4111.

**[9]** Shakya, M.; Varum, H.; Vicente, R.; Costa, A. Empirical formulation for estimating the fundamental frequency of slender masonry structures. International Journal of Architectural Heritage 2016, 10, 55–66. https://doi.org/https://doi.org/10.1080/15583058.2014.951796.

**[10]** Bartoli, G.; Betti, M.; Marra, A.M.; Monchetti, S. On the role played by the openings on the first frequency of historic masonry towers. Bulletin of Earthquake Engineering 2020, 18, 427–451. https://doi.org/https://doi.org/10.1007/s10518-019-00662-9.

**[11]** Pallarés, F.J.; Betti, M.; Bartoli, G.; Pallarés, L. Structural health monitoring (SHM) and Nondestructive testing (NDT) of slender masonry structures: A practical review. Construction and Building Materials 2021, 297, 123768. https://doi.org/https://doi.org/10.1016/j.conbuildmat.2021.123768.
