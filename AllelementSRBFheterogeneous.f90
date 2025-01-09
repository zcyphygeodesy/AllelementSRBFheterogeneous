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
      !����������ɢ�۲ⳡԪ�ļ���¼����Լ����ʽ�����/վ�������ȣ���С������γ�ȣ���С��������ظߣ�m�����в�۲ⳡԪ������
        !��Ԫ���ͣ�0~5����Ȩֵ��������¼ǰ5�����Ե�λ�ú�˳��Լ�����䡣
      !The agreed format of the observation file record: ID (point no / station name), longitude (degree decimal), 
       !latitude, ellipsoidal height (m), observation, ..., observation type (0 ~ 5), weight, ... The order of the 
       !first five attributes is fixed by convention.
	write(observationfl,*) "obsrgagrr.txt" !The discrete heterogeneous observation file
	write(surfhgtgrdfl,*) "surfhgt.dat"!The ellipsoidal height grid file of the calculation surface
      para(1)=360;para(2)=1800!SRBF��С������ Minimum and maximum degree
      para(3)=3!�༶���� the order number m
      para(4)=0!0-����༶�Ӻ˺�����1-PoissonС���˺���
      para(5)=1800!SRBF��Reuter�ȼ�K the Reuter network level K for the SRBFs
      para(6)=100.d0!SRBF�������þ���(km) the action distance of SRBF center
      para(7)=10.d0!Bjerhammar������Ե����ƽ����� the Bjerhammar sphere burial depth (km)
      !�۲ⳡԪ������Ȩֵ����� the column ordinal number of the type and weight of the observation
      !=0�Ŷ�������=1�߳��쳣��=2�ռ��쳣��=3�Ŷ������ݶȣ�=4����ƫ��
      !=0 gravity disturbance (mGal), =1 height anomaly (m), =2 gravity anomaly (mGal), =3 disturbing gravity gradient
        !(E, radial) or =4 vertical deflection (��)
      para(8)=6;para(9)=7
      !��Ȩֵ���������С��1���򳬳���¼����ţ����ļ���¼��Ȩֵ����С����ʱ������Ĭ�ϵ�Ȩ��
      !���ļ���¼��Ȩֵ������ʱ���ù۲���������SRBFϵ�����ƣ��������н����ֱ�����ڲⶨ�ù۲������ⲿ����ָ�ꡣ
      !When the column ordinal number of the weight attribute is less than 1, exceeds the column number of the record, 
        !or the weight is less than zero, the program makes the weight equal to 1.
      !When the weight in the file record is equal to zero, the observation will not participate in the estimation of the 
        !SRBF coefficients, and the program can be employed to measure the external accuracy index of the observations.
      !ѡ��ɵ��س�Ԫ���ͣ����ÿɵ��ع۲ⳡԪ�����ʦʡ�
      !Select the type of the adjustable observations and set the contribution rate �� of the adjustable observations.
      para(10)=0;para(11)=1.d0
      !ѡ�񷨷��̽��㷽�� Select the method of the solution of normal equation
      para(12)=1!1-LU�ֽ�,2-Cholesky�ֽ�,3-��С����QR�ֽ�,4-��С��������ֵ�ֽ�,5-�����
      write(*, *)"    Begin compulation, please wait......"
      call SRBFheterogeneous(observationfl,surfhgtgrdfl,para)
      write (*,*)'    Complete the computation! All approached target field elements are saved'
      write (*,*)'      in the file SRBFhetero.txt.'
      write (*,*)'    The program outputs also the residual observation file residuals.txt into'
      write (*,*)'      the current directory.'
      pause
      end
