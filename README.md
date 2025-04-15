![logo](https://github.com/CMCC-Foundation/OceanVar2/blob/main/doc/logo.png)

# OceanVar2 

**OceanVar2,** developed at CMCC, is an open-source state-of-the-art variational ocean data assimilation framework. [![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
**OceanVar2** is used to study the ocean across a range of scales, from global to local and short-term to multi-decadal. Additionally, **OceanVar2** provides initial conditions for short-range to seasonal forecasts.

**OceanVar** was first introduced by [Dobricic and Pinardi (2008)](https://www.sincem.unibo.it/images/articoli/10.1016_j.ocemod.2008.01.004.pdf)  and is based on a three-dimensional variational method, formulated in its classical incremental variant. **OceanVar2** features a modular design that allows for flexibility in incorporating diverse data sources and error covariance representations. It supports the assimilation of in-situ ocean measurements of temperature and salinity, as well as remotely sensed data, including altimetry, and sea surface temperature and salinity, from a variety of datasets depending on the specific application.

**OceanVar2** decomposes the background error covariance matrix into physically based linear operators, allowing for individual analysis of specific error components. A key feature of **OceanVar2** is its ability to represent error correlations between temperature, salinity, and sea level anomalies (the sea level operator). OceanVar2 offers the flexibility of using either a dynamic height or a barotropic model for closed domains, or EOF-based correlations. Horizontal error correlations are modeled with a diffusive operator, replacing the previous recursive filter. Furthermore, the **OceanVar2** code has been extensively revised into a unified, consistent, and fully parallelized framework, integrating past developments.
##

# [user guide](https://github.com/CMCC-Foundation/OceanVar2/blob/main/doc/OceanVar_User_Manual.pdf)
PDF User guide can be downloaded [here](https://github.com/CMCC-Foundation/OceanVar2/blob/main/doc/OceanVar_User_Manual.pdf).
Doxigen documentation can be generated using the [Doxyfile](https://github.com/CMCC-Foundation/OceanVar2/blob/main/doc/Doxyfile)

# [Test Case](https://github.com/CMCC-Foundation/MedFS831)

To help users become familiar with the code, a test-case is also provided [here](https://github.com/CMCC-Foundation/MedFS831)


## References using OceanVar:

 - **2006**:
    *  [Dobricic, S., Pinardi, N., Adani, M., Bonazzi, A., Fratianni, C. and Tonani, M.Mediterranean Forecasting System: An improved assimilation scheme for sea-level anomaly and its validation. Q. J. R. Meteorol. Soc. 131: 3627–3642](https://www.sincem.unibo.it/images/articoli/10.1256_qj.05.100.pdf) 

 - **2007**:
    *  [Dobricic, S., Pinardi, N., Adani, M., Tonani, M., Fratianni, C., Bonazzi, A., Fernandez, V. Daily oceanographic analyses by the Mediterranean basin scale assimilation system. Ocean
   Sciences 3, 149–157](https://os.copernicus.org/articles/3/149/2007/os-3-149-2007.pdf)
   
-  **2008**:
    *   [Dobricic, S. and Pinardi, N.: An oceanographic three-dimensional variational data assimilation scheme, Ocean Model.,    22, 89–105, 2008.](https://www.sincem.unibo.it/images/articoli/10.1016_j.ocemod.2008.01.004.pdf) 
   
 - **2011**
    * [Adani, M., S. Dobricic, and N. Pinardi: Quality Assessment of a 1985–2007 Mediterranean Sea Reanalysis. J. Atmos. Oceanic    Technol., 28, 569–589](https://doi.org/10.1175/2010JTECHO798.1)
    * [Storto, A., S. Dobricic, S. Masina, and P. Di Pietro: Assimilating Along-Track Altimetric Observations through Local    Hydrostatic Adjustment in a Global Ocean Variational Assimilation    System. Mon. Wea. Rev., 139, 738–754](https://doi.org/10.1175/2010MWR3350.1).
   
- **2012**:
    * [Dobricic, S., Dufau, C.,    Oddo, P., Pinardi, N., Pujol, I., and Rio, M.-H.: Assimilation of SLA    along track observations in the Mediterranean with an oceanographic    model forced by atmospheric pressure, Ocean Sci., 8, 787–795](https://doi.org/10.5194/os-8-787-2012)
    * [Nilsson, J. A. U.,    Dobricic, S., Pinardi, N., Poulain, P.-M., and Pettenuzzo, D.
   Variational assimilation of Lagrangian trajectories in the    Mediterranean ocean Forecasting System Ocean Sci., 8, 249-259](http://dx.doi.org/10.5194/os-8-249-2012)
   
- **2014**:
    * [Teruzzi, A., S. Dobricic,    Solidoro, C., Cossarini, G. A 3-D variational assimilation scheme in    coupled transport-biogeochemical models: forecast of Mediterranean    biogeochemical properties. J. Geophys. Res.: Oceans, 119 (1)](https://doi.org/10.1002/2013JC009277).
   
- **2015**:
    * [Dobricic, S., Wikle, C.    K., Milliff, R. F., Pinardi, N., Berliner, L. M. Assimilation of    oceanographic observations with estimates of vertical    background-error covariances by a Bayesian hierarchical model Quarterly Journal of the Royal Meteorological Society, 141, 182-194](http://dx.doi.org/10.1002/qj.2348)
   
- **2016**:
    * [Oddo, P. Storto, A. Dobricic, S. Russo, A. Lewis, C. Onken, R. Coelho, E., A hybridvariational-ensemble data assimilation scheme with systematic error correction for limited-area ocean models, Ocean Science,12,2016,5,1137-1153](https://os.copernicus.org/articles/12/1137/2016/os-12-1137-2016.pdf)
    * [Aydoğdu, A., Pinardi, N.,    Pistoia, J., Martinelli, M., Belardinelli, A., Sparnocchia, S.,Assimilation experiments for the Fishery Observing System in the    Adriatic Sea, Journal of Marine Systems, Volume 162, October 2016,    Pages 126-136, ISSN 0924-7963, ](https://doi.org/10.1016/j.jmarsys.2016.03.002 )
    * [Storto, A., Masina, S. and    Navarra, A., Evaluation of the CMCC eddy-permitting global ocean    physical reanalysis system (C-GLORS, 1982–2012) and its assimilation components. Q.J.R. Meteorol. Soc., 142: 738-758.](https://doi.org/10.1002/qj.2673)
    * [Storto, A., Variational    quality control of hydrographic profile data with non-Gaussian errors    for global ocean variational data assimilation systems, Ocean    Modelling, Volume 104, 2016, Pages 226-241, ISSN 1463-5003](https://doi.org/10.1016/j.ocemod.2016.06.011).
   
- **2018**:
    * Storto, A., Oddo, P.,    Cipollone, A., Mirouze, I., Lemieux-Dudon, B. Extending an    oceanographic variational scheme to allow for affordable hybrid and    four-dimensional data assimilation. Ocean Modelling, 128, pp. 67-86, [10.1016/j.ocemod.2018.06.005](https://doi.org/10.1016/j.ocemod.2018.06.005)
    * Teruzzi, A., Bolzon, G.,    Salon, S., Lazzari, P., Solidoro, C., Cossarini, G. Assimilation of    coastal and open sea biogeochemical data to improve phytoplankton    simulation in the mediterranean sea. Ocean Model., 132 (2018), pp.    46-60
   
- **2019**:
    * Storto A, Oddo P. Optimal    Assimilation of Daytime SST Retrievals from SEVIRI in a Regional    Ocean Prediction System. Remote Sensing. 2019; 11(23):2776.https://doi.org/10.3390/rs11232776
- **2020**:
    * Cipollone A., A. Storto &    S. Masina, "Implementing a parallel version of a variational scheme    in a global assimilation system at eddy-resolving resolution",JTECH,37(10),1865-1876,    [https://doi.org/10.1175/JTECH-D-19-0099.1](https://doi.org/10.1175/JTECH-D-19-0099.1)
    * Storto, A., Falchetti, S.,    Oddo, P., Jiang, Y.-M., & Tesei, A. Assessing the impact of different    ocean analysis schemes on oceanic and underwater acoustic    predictions. Journal of Geophysical Research: Oceans, 125,    e2019JC015636.
   [https://doi.org/10.1029/2019JC015636](https://doi.org/10.1029/2019JC015636)
- **2021**:
    * Storto, A., G. De    Magistris, S. Falchetti, and P. Oddo: A Neural Network–Based    Observation Operator for Coupled Ocean–Acoustic Variational Data    Assimilation. Mon. Wea. Rev., 149, 1967–1985,
   [https://doi.org/10.1175/MWR-D-20-0320.1](https://doi.org/10.1175/MWR-D-20-0320.1)
    * Escudier, R., Clementi,    E., Cipollone, A., Pistoia, J., Drudi, M., Grandi, A., Lyubartsev,    V., Lecci, R., Aydogdu, A., Delrosso, D., et al.: A high-resolution    reanalysis for the Mediterranean Sea, Frontiers in Earth Science, 9,    702285, 2021.
    * Lima L., Ciliberti S.A.,    Aydogdu A., Masina S., Escudier R., Cipollone A., Azevedo D., Causio    S., Peneva E., Lecci R., Clementi E., Jansen E., licak M., Creti S.,    Stefanizzi L., Palermo F. & Coppini G., "Climate signals in the black    sea from a multidecadal eddy- resolving Reanalysis",Front. Marine    Sci. 8, 1214 , [https://doi.org/10.3389/fmars.2021.710973](https://doi.org/10.3389/fmars.2021.710973)
- **2022**:
    * Ciliberti, S.A.; Jansen,    E.; Coppini, G.; Peneva, E.; Azevedo, D.; Causio, S.; Stefanizzi, L.;    Creti, S.; Lecci, R.; Lima, L.; et al. The Black Sea Physics Analysis    and Forecasting System within the Framework of the Copernicus Marine    Service. J. Mar. Sci. Eng. 2022, 10, 48.    [https://doi.org/10.3390/jmse10010048](https://doi.org/10.3390/jmse10010048)
    * Oddo P, Falchetti S, Viola    S, Pennucci G, Storto A, Borrione I, Giorli G, Cozzani E, Russo A,    Tollefsen C. Evaluation of different Maritime rapid environmental    assessment procedures with a focus on acoustic performance. J Acoust    Soc Am. 2022 Nov;152(5):2962. doi: 10.1121/10.0014805. PMID: 36456253.
    * Cipollone A., Banerjee, D.    S., Iovino, D., Aydogdu, A., and Masina, S. (2023): "Bivariate sea-ice assimilation for global-ocean analysis-reanalysis", Ocean    Sci.,19,1375-1392, https://doi.org/10.5194/egusphere-2022-1337
- **2023**:
  * Clementi, E., Drudi, M.,    Aydogdu, A., Moulin, A., Grandi, A., Mariani, A., Goglio, A. C.,    Pistoia, J., Miraglio, P., Lecci, R., Palermo, F., Coppini, G.,    Masina, S., & Pinardi, N. (2023). Mediterranean Sea Physical Analysis    and Forecast (CMS MED-Physics, EAS8 system) (Version 1) [Data set]. Copernicus Marine Service (CMS).[https://doi.org/10.25423/CMCC/MEDSEA_ANALYSISFORECAST_PHY_006_013_EAS8](https://doi.org/10.25423/CMCC/MEDSEA_ANALYSISFORECAST_PHY_006_013_EAS8)
  * Coppini, G., Clementi, E.,    Cossarini, G., Salon, S., Korres, G., Ravdas, M., Lecci, R., Pistoia,    J., Goglio, A. C., Drudi, M., et al.: The Mediterranean forecasting    system. Part I: evolution and performance, EGUsphere, pp. 1–50, 2023.
