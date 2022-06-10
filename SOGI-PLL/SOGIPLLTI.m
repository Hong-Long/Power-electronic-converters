clear all;
close all;
clc;
% define the math type being used on the controller using objects from the
% fixed-pointtoolbox in matlab
% Select numeric type,let's choose Q23
T = numerictype('WordLength',32,'FractionLength',23);
% Specify mathattributes to the fimath object
F=fimath('RoundMode','floor','OverflowMode','wrap');
F.ProductMode='SpecifyPrecision';
F.ProductWordLength=32;
F.ProductFractionLength=23;
F.SumMode='SpecifyPrecision';
F.SumWordLength=32;
F.SumFractionLength=23;
%specify fipref object,to display warning in cases of overflowand
%under flow
P=fipref;
P.LoggingMode='on';
P.NumericTypeDisplay='none';
P.FimathDisplay='none';
%PLL Modelling starts from here
Fs=50000;%Sampling frequency= 50Khz
GridFreq=50;%Nominal Grid Frequency in Hz
Tfinal=0.2;%Time the simulationis run for = 0.5 seconds
Ts=1/Fs;%Sampling Time= 1/Fs
t=0:Ts:Tfinal;%Simulation Time vector
wn=2*pi*GridFreq;%Nominal Grid Frequency in radians
%declare arrays used by the PLL process
err=fi([0,0,0,0,0],T,F);
ylf=fi([0,0,0,0,0],T,F);
Mysin=fi([0,0,0,0,0],T,F);
Mycos=fi([1,1,1,1,1],T,F);
theta=fi([0,0,0,0,0],T,F);
dc_err=fi([0,0,0,0,0],T,F);
wo=fi(0,T,F);

% used for plotting
Plot_Var=fi([0,0,0,0],T,F);
Plot_theta=fi([0,0,0,0],T,F);
Plot_osgu=fi([0,0,0,0],T,F);
Plot_osgqu=fi([0,0,0,0],T,F);
Plot_D=fi([0,0,0,0],T,F);
Plot_Q=fi([0,0,0,0],T,F);
Plot_dc_err=fi([0,0,0,0,0],T,F);
% orthogonal signal generator
% using trapezoidal approximation
osg_k=0.5;
osg_x=2*osg_k*wn*Ts;
osg_y=(wn*wn*Ts*Ts);
osg_b0=osg_x/(osg_x+osg_y+4);
osg_b2=-1*osg_b0;
osg_a1=(2*(4-osg_y))/(osg_x+osg_y+4);
osg_a2=(osg_x-osg_y-4)/(osg_x+osg_y+4);

osg_qb0=(osg_k*osg_y)/(osg_x+osg_y+4);
osg_qb1=2*osg_qb0;
osg_qb2=osg_qb0;

osg_k=fi(osg_k,T,F);
osg_x=fi(osg_x,T,F);
osg_y=fi(osg_y,T,F);
osg_b0=fi(osg_b0,T,F);
osg_b2=fi(osg_b2,T,F);
osg_a1=fi(osg_a1,T,F);
osg_a2=fi(osg_a2,T,F);
osg_qb0=fi(osg_qb0,T,F);
osg_qb1=fi(osg_qb1,T,F);
osg_qb2=fi(osg_qb2,T,F);
osg_u=fi([0,0,0,0,0,0],T,F);
osg_qu=fi([0,0,0,0,0,0],T,F);

u_Q=fi([0,0,0],T,F);
u_D=fi([0,0,0],T,F);
% generate input signal

% CASE1 : Phase Jump at the Mid Point
L=length(t);
for n=1:floor(L)
    u(n)=sin(2*pi*GridFreq*Ts*n);
end
for n=1:floor(L)
    u1(n)=sin(2*pi*GridFreq*Ts*n);
end
for n=floor(L/2):L
    u(n)=sin(2*pi*GridFreq*Ts*n+pi/2);
end
u=fi(u,T,F);
% simulate the PLL process 
for n=3:Tfinal/Ts % No of iteration of the PLL processin the simulation time
    %OrthogonalSignalGenerator
    osg_u(1)=(osg_b0*(u(n)-u(n-2)))+osg_a1*osg_u(2)+osg_a2*osg_u(3);
    osg_u(3)=osg_u(2);
    osg_u(2)=osg_u(1);
    osg_qu(1)=(osg_qb0*u(n)+osg_qb1*u(n-1)+osg_qb2*u(n-2))+osg_a1*osg_qu(2)+osg_a2*osg_qu(3);
    osg_qu(3)=osg_qu(2);
    osg_qu(2)=osg_qu(1);
    % park trasnform from alpha beta to d-q axis
    u_Q(1)=Mycos(2)*osg_u(1)+Mysin(2)*osg_qu(1);
    u_D(1)=-Mysin(2)*osg_u(1)+Mycos(2)*osg_qu(1);
    
    %Loop Filter
    ylf(1)=fi(1,T,F)*ylf(2)+fi(166.877556,T,F)*u_Q(1)+fi(-166.322444,T,F)*u_Q(2);
    
    u_Q(2)=u_Q(1);
    u_D(2)=u_D(1);
    %Limit LF accordingto its Q? sizepipeline
    ylf(1)=max([ylf(1) fi(-128,T,F)]);
    ylf(1)=min([ylf(1) fi(128,T,F)]);
    ylf(2)=ylf(1);
    
    %update output frequency
    wo=GridFreq+ylf(1);
    
    %update the output phase
    theta(1)=theta(2)+wo*fi(Ts,T,F);
    if(theta(1)>fi(1.0,T,F))
        theta(1)=fi(0,T,F);
    end
    theta(2)=theta(1);
    Mysin(1)=sin(theta(1)*fi(2*pi,T,F));
    Mycos(1)=cos(theta(1)*fi(2*pi,T,F));
    Mysin(2)=Mysin(1);
    Mycos(2)=Mycos(1);
    
    Plot_theta(n+1)=theta(1);
    Plot_osgu(n+1)=osg_u(1);
    Plot_osgqu(n+1)=osg_qu(1);
    Plot_Var(n+1)=Mysin(1);
    Plot_D(n+1)=u_D(1);
    Plot_Q(n+1)=u_Q(1);
end
% CASE1 : Phase Jump at the Mid Point
% error=Plot_Var-u;

%CASE2 : Harmonics
%error=Plot_Var-u1;

%CASE3: FrequencyVariations
error=Plot_Var-u;

%CASE4: AmplitudeVariations
%error=Plot_Var-u1;
subplot(3,1,1),plot(t,Plot_Var,'r',t,u,'b'),title('SPLL(red)& IdealGrid(blue)');
subplot(3,1,2),plot(t,error,'r'),title('Error');
subplot(3,1,3),plot(t,u1,'r',t,Plot_Var,'b'),title('SPLLOut(Blue)& IdealGrid(Red)');
