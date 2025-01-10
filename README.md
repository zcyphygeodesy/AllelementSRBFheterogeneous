## Fortran codes for all-element field modelling using SRBFs from heterogeneous observations
https://www.zcyphygeodesy.com/en/h-nd-155.html
### [Algorithm purpose]
&emsp;```From various heterogeneous observations which can be the residual gravity disturbance, height anomaly, gravity anomaly, disturbing gravity gradient, or vertical deflection, determinate the residual gravity disturbance, height anomaly, gravity anomaly, disturbing gravity gradient, and vertical deflection outside geoid using spherical radial basis functions (SRBFs), to realize the unified modelling on regional gravity field and geoid.```  
&emsp;```From various heterogeneous observations, determinate the all-element residual gravity field models outside geoid using SRBFs, to realize the unified modelling on regional gravity field and geoid.```  
&emsp;```The program is a high performance and adaptable modelling tool on local gravity field. Various observations with heterogeneity, different altitudes, cross-distribution and land-sea coexisting can be directly employed to estimate the all-element models of gravity field without reduction, continuation and gridding.```  
&emsp;```The program has strong capacity on the detection of observation gross errors, measurement of external accuracy indexes and control of computation performance.```
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9OzltwYo0PbRlAQwpQ047gg.jpg)
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9ezltwYo6LHlwQYwpQ047gg.jpg)
### [Main program for test entrance]
    AllelementSRBFheterogeneous.f90
    Input parameter: observationfl - The heterogeneous observation file name.
      ▪ The agreed format of the observation file record: ID (point no / station name), longitude (degree decimal), latitude, ellipsoidal height (m), observation, ..., observation type (0 ~ 5), weight, ... The order of the first five attributes is fixed by convention.
      ▪ The observation types and units: 0 - residual gravity disturbance (mGal), 1 - residual height anomaly (residual geoidal height, m), 2 - residual gravity anomaly (mGal), 3 - residual disturbing gravity gradient (E, radial), 4 - residual vertical deflection southward component (″), 5 - residual vertical deflection westward component (″).
    Input parameter: surfhgtgrdfl - The ellipsoidal height grid file name of the calculation surface. There is no limit to the grid resolution of the calculation surface.
    Input parameters: para(1:7) - the minimum and maximum degree of SRBF Legendre expansion, the order number m, the spherical radial basis functions (=0 radial multipole kernel function, =1 Poisson wavelet kernel function), Reuter network level K, action distance (km) of SRBF center and Bjerhammar sphere burial depth (km).
    Input parameters: para(8:9) - the column ordinal numbers of the types of the observation and weight.
      ▪ When the column ordinal number of the weight attribute is less than 1, exceeds the column number of the record, or the weight is less than zero, the program makes the weight equal to 1.
      ▪ When the weight in the file record is equal to zero, the observation will not participate in the estimation of the SRBF coefficient, and the program can be employed to measure the external accuracy index of the observations.
    Input parameters: para(10:11) - the type of the adjustable observations and the contribution rate of the adjustable observations.
    Input parameters: para(12) - method of the solution of normal equation. =1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square QR decomposition, =4 Minimum norm singular value decomposition and =5 ridge dstimation.
### (1) Module for all-element modelling on gravity field using SRBFs
    SRBFheterogeneous(observationfl,surfhgtgrdfl,para)
    Input parameter: observationfl - The heterogeneous residual anomalous field element observation file name.
    Input parameter: surfhgtgrdfl - The ellipsoidal height grid file name of the calculation surface.
    Output the all-element approached residual field element file SRBFhetero.txt.
      ▪ The file record: ID, longitude (degree decimal), latitude, ellipsoidal height (m) of calculated point, residual gravity disturbance (mGal), residual height anomaly (m), residual gravity anomaly (mGal), residual disturbing gravity gradient (E, radial) and residual vertical deflection (″, SW).
The module outputs the all-element grid files into the current directory, whose grid specification are the same as the input ellipsoidal height grid of calculation surface.
      ▪ The residual gravity disturbance grid file SRBFhetero.rga (mGal)
      ▪ residual height anomaly grid file SRBFhetero.ksi (m)
      ▪ residual gravity anomaly grid file SRBFhetero.gra (mGal)
      ▪ residual disturbing gravity gradient grid file SRBFhetero.grr (E, radial) and 
      ▪ residual vertical deflection vector grid file SRBFhetero.dft (″, SW). 
    The module also outputs the residual observation file residuals.txt.
      ▪ The statistical results of each type of observations occupies a row of file header, whose format: observation type (0~5), source observation mean, standard deviation, minimum, maximum, residual observation mean, standard deviation, minimum, maximum. 
      ▪ The record format: ID, longitude, latitude, ellipsoidal height, residual observation, source observation, observation type, weight..
The module also outputs the SRBF spatial curve file SRBFspc.txt, spectral curve file SRBFdgr.txt of 5 kinds of field elements and SRBF center file SRBFcenter.txt into the current directory.
      ▪ SRBFspc.txt file header format: SRBF type (0-radial multipole kernel function, 1-Poisson wavelet kernel function), order of SRBF, Minimum and maximum degree of SRBF Legendre expansion, Bjerhammar sphere buried depth (km). The record format: spherical distance (km), the normalized SRBF values from the gravity disturbance, height anomaly, disturbing gravity gradient and total vertical deflection.
      ▪ The file header of SRBFdgr.txt is the same as SRBFspc.txt. The record format: the degree n of SRBF Legendre expansion, degree n normalized SRBF values from the gravity disturbance, height anomaly, disturbing gravity gradient and total vertical deflection.
      ▪ SRBFcenter.txt file header format: Reuter grid level, SRBF center number, cell-grid number in meridian circle direction, maximum cell-grid number in prime vertical circle direction, latitude interval ('). The record format: point no, longitude (degree decimal), geocentric latitude, cell-grid area deviation percentage, longitude interval of cell-grid in prime vertical circle direction (').
### (2) Computation module for the Reuter network parameters
    ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)
    Input parameters: rhd(4) - minimum and maximum longitude, minimum and maximum geocentric latitude of the Reuter network.
    Input parameters: lvl, nn, mm - the Reuter network level, maximum number of Reuter centers in the meridian direction and that in the parallel direction.
    Return parameter: Blat - the geocentric latitude (degree decimal) of Reuter centers in the first parallel direction.
    Return parameters: Kt - the number of Reuter centers, equal to the number of unknowns to be estimated.
    Return parameters: nln(nn) - the number of the Reuter centers in the parallel direction.
    Return parameters: sr(nn) - the percentage of the difference between the area of the Reuter cell-grid in the parallel direction and the area of the equatorial cell-grid.
    Return parameters: dl(nn) - the longitude interval (degree decimal) of the Reuter centers in the parallel direction.
    Return parameters: nrd(nn,mm) - ordinal number value of the Reuter centers.
    Return parameters: lon(nn,mm) - longitude value of the Reuter centers.
### (3) Module for the best match between the Reuter centers and observation points
    Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
    The module calculates the number of observation points in the Reuter cell-grid and the number gpnt of Reuter centers corrected.
    Return parameter: edgn - the number of Reuter centers around the edge of the Reuter grid.
    Return parameters: enode(edgn) - the ordinal number value of Reuter centers in the edge of the Reuter grid.
    Return parameter: rlatlon(edgn,2) - the geocentric latitude and longitude (degree decimal) of Reuter centers in the edge of the Reuter grid.
### (4) Computation module for the SRBF curves of 5 kinds of field elements
    SRBF5all(RBF,order,krbf,mpn,mdp,minN,maxN,NF,nta)
    Return parameters: RBF(NF+1,5) - The SRBF curves of 5 kinds of field elements, which are calculated by the action distance of SRBF center and Reuter grid level.
      ▪ Where RBF(NF+1,knd): knd=1 gravity disturbance (mGal), =2 height anomaly (m), =3 gravity anomaly (mGal), =4 disturbing gravity gradient (E, radial) or =5 vertical deflection (″).
    Input parameters: nta - the bandwidth parameter and nta = (r0-dpth)/r0, here dpth is the Bjerhammar sphere burial depth and r0 is the average geocentric distance of observation points.
    Input parameters: krbf,order - krbf=0 radial multipole kernel function, =1 Poisson wavelet kernel function and order is the order number of SRBF.
    Input parameters: mpn(maxN-minN+1, NF+1), mdp(maxN-minN+1, NF+1) - all minN to maxN degree Legendre functions and their first derivatives.
### (5) Calculation module for the position of the calculation point in the Reuter grid
    RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)
    Input parameters: rln(3) - the spherical coordinates of the calculation point。
    Return parameters: ki,kj - the position of the calculation point rln(3) in the Reuter grid, It is represented by the element of the 2-D ordinal number array of the Reuter grid and ki>0, kj>0.
### (6) Solution module of large positive definite symmetric equations
    Equsolve(BB,xx,nn,BL,knd,bf)
    The mkl_lapack95_ilp64.lib library is called to solve the large equations BB.xx = BL. bf(8) is the property of the solution.
    Input parameter: knd - method of the solution of normal equation and knd =1 LU triangular decomposition method, =2 Cholesky decomposition, =3 least square QR decomposition, =4 Minimum norm singular value decomposition.
### (7) Calculation module for the normal gravity field
    normdjn(GRS,djn); GNormalfd(BLH,NFD,GRS)
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass).
### (8) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt1(pn,dp1,n,t) ! t=cos ψ
### (9) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); RLAT_BLH(GRS, RLAT, BLH)
### (10) Other auxiliary modules
    LegPn02(mpn,mdp,mp2,minN,maxN,NF,dr); PickRecord(line, kln, rec, nn)
    RBFvalue(RBF(:,1),NF,dr,dln(2),tmp); drln(rln,rlnk,dln); Stat1d(dt,nn,rst)
### [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler. mkl_lapack95_ilp64.lib link library required.
### [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    7.10 Theory and algorithm of gravity field approach using spherical radial basis functions
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file and all input and output data.
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9OzltwYopYzL-AYwpQ047gg.jpg)
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9ezltwYoiP3aaDClDTjuCA.jpg)
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9uzltwYouLrdjQQwpQ047gg.jpg)
![](https://24192633.s21i.faiusr.com/2/ABUIABACGAAg9uzltwYozIPxeTClDTjuCA.jpg)
