clear
clc

clockstart = clock

M=1.7565670*10.^28;%%Mass of sun
v= 0; %%v=true anomaly (angle from periapsis to planet)
%%v=2*atan(sqrt((1+e)/1-e)*tan(E));
G=6.67384*10.^-11;
T=zeros(0,2000*7+2);%%Time
R=zeros(0,2000*7+2);%%Radius
ang=zeros(0,2000*7+2);%%angle
h = zeros(0,2000*7+2);

i=zeros(0,2000*7+2);
count_1=-1;
%%semi-major axis a
for planet=1:1:7
t=0;

 A=[1, 5263138304; 2, 9832684544; 3, 13599840256; 4, 20726155264; 5, 40839348203; 6, 68773560320; 7, 90118820000];
 e_e=[1, .2; 2, .01; 3, 0; 4, .05; 5, .14; 6, .05; 7, .26];
 I_i=[1,7; 2, 2.1; 3, 0; 4, .06; 5, 5; 6, 1.304; 7, 6.15];
 m_0=[1, pi; 2, pi; 3, pi; 4, pi; 5, pi; 6, .1; 7, pi]
% planet = input('what planet?     ');
a=A(planet,2);
e=e_e(planet,2);
I=I_i(planet,2);
%   a=input('What is the semi-major axis?     ');%%Semi-Major Axis
%   e=input('what is the eccentricity?     ');%%Eccentricity
%   I=input('what is the inclination?     ');%%Incliination


count_1=count_1+1;
count=0;
%%i = tan(I+(pi/2));
I=I.*(2*pi/360);


for t=0:156992048/2000:156992048
syms E;
n=sqrt((G*M)/(a.^3));
count=count+1
    T(count)=t;
m=mod(n*(t)+m_0(planet,2), 2.*pi);
E=solve(m==E-e*sin(E));
%%v=2*atan(sqrt(1-e)*cos(E/2)/(sqrt(1+e)*sin(E/2)));
%%v=acos((cos(E-e))/(1+e*cos(E)));
v=mod(double(2*atan(sqrt((1+e)/(1-e))*tan(E/2))+2*pi) ,2*pi);

    r=(a*(1-e.^2))/(1+e.*cos(v));
    R(500*count_1+count)=r;
    ang(500*count_1+count)=v;
    i_mod=-I*sin(v+pi/2);
    i(500*count_1+count)=i_mod;
    H=r*sin(i_mod);
    h(500*count_1+count)=H;

end
end
    consolidated = [T R ang];
[X,Y,Z]=pol2cart(ang,R,h);
plot3(X,Y,Z)
%%plot(T,ang)
endtime = clock
time_elapsed = etime(endtime,clockstart)
