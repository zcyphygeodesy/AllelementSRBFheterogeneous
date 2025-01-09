!  AllelementSRBFheterogeneous.f90 
!
!  FUNCTIONS:
!  AllelementSRBFheterogeneous - Entry point of console application.
!
!****************************************************************************
      program AllelementSRBFheterogeneous
      implicit none
	character*800::observationfl,surfhgtgrdfl
      real*8 para(20)
!---------------------------------------------------------------------
      !多种异质离散观测场元文件记录采用约定格式：点号/站名，经度（度小数），纬度（度小数），大地高（m），残差观测场元，…，
        !场元类型（0~5），权值，…。记录前5项属性的位置和顺序约定不变。
      !The agreed format of the observation file record: ID (point no / station name), longitude (degree decimal), 
       !latitude, ellipsoidal height (m), observation, ..., observation type (0 ~ 5), weight, ... The order of the 
       !first five attributes is fixed by convention.
	write(observationfl,*) "obsrgagrr.txt" !The discrete heterogeneous observation file
	write(surfhgtgrdfl,*) "surfhgt.dat"!The ellipsoidal height grid file of the calculation surface
      para(1)=360;para(2)=1800!SRBF最小最大阶数 Minimum and maximum degree
      para(3)=3!多级次数 the order number m
      para(4)=0!0-径向多级子核函数，1-Poisson小波核函数
      para(5)=1800!SRBF的Reuter等级K the Reuter network level K for the SRBFs
      para(6)=100.d0!SRBF中心作用距离(km) the action distance of SRBF center
      para(7)=10.d0!Bjerhammar球面相对地面的平均深度 the Bjerhammar sphere burial depth (km)
      !观测场元类型与权值列序号 the column ordinal number of the type and weight of the observation
      !=0扰动重力，=1高程异常，=2空间异常，=3扰动重力梯度，=4垂线偏差
      !=0 gravity disturbance (mGal), =1 height anomaly (m), =2 gravity anomaly (mGal), =3 disturbing gravity gradient
        !(E, radial) or =4 vertical deflection (″)
      para(8)=6;para(9)=7
      !当权值属性列序号小于1，或超出记录列序号，或文件记录中权值属性小于零时，程序默认等权。
      !当文件记录中权值等于零时，该观测量不参与SRBF系数估计，程序运行结果可直接用于测定该观测量的外部精度指标。
      !When the column ordinal number of the weight attribute is less than 1, exceeds the column number of the record, 
        !or the weight is less than zero, the program makes the weight equal to 1.
      !When the weight in the file record is equal to zero, the observation will not participate in the estimation of the 
        !SRBF coefficients, and the program can be employed to measure the external accuracy index of the observations.
      !选择可调控场元类型，设置可调控观测场元贡献率κ。
      !Select the type of the adjustable observations and set the contribution rate κ of the adjustable observations.
      para(10)=0;para(11)=1.d0
      !选择法方程解算方法 Select the method of the solution of normal equation
      para(12)=1!1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4-最小范数奇异值分解,5-岭估计
      write(*, *)"    Begin compulation, please wait......"
      call SRBFheterogeneous(observationfl,surfhgtgrdfl,para)
      write (*,*)'    Complete the computation! All approached target field elements are saved'
      write (*,*)'      in the file SRBFhetero.txt.'
      write (*,*)'    The program outputs also the residual observation file residuals.txt into'
      write (*,*)'      the current directory.'
      pause
      end
