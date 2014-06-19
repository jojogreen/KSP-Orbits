clear
clc

clockstart = clock

M=1.7565670*10.^28;%%Mass of sun
v= 0; %%v=true anomaly (angle from periapsis to planet)
%%v=2*atan(sqrt((1+e)/1-e)*tan(E));
G=6.67384*10.^-11;
%%Time
num_sec= 6000;
T=zeros(0,num_sec+3);
%%Radius for planets
R_1=zeros(0,num_sec+3);%%Radius
R_2=zeros(0,num_sec+3);
R_3=zeros(0,num_sec+3);
R_4=zeros(0,num_sec+3);
R_5=zeros(0,num_sec+3);
R_6=zeros(0,num_sec+3);
R_7=zeros(0,num_sec+3);

%%angle for planets
ang=zeros(0,num_sec+3);%%angle
ang=zeros(0,num_sec+3);
ang=zeros(0,num_sec+3);
ang=zeros(0,num_sec+3);
ang=zeros(0,num_sec+3);
ang=zeros(0,num_sec+3);
ang=zeros(0,num_sec+3);

%%Heights for planets
h_1= zeros(0,num_sec+3);
h_2= zeros(0,num_sec+3);
h_3= zeros(0,num_sec+3);
h_4= zeros(0,num_sec+3);
h_5= zeros(0,num_sec+3);
h_6 = zeros(0,num_sec+3);
h_7 = zeros(0,num_sec+3);

i_1=zeros(0,num_sec+3);
i_2=zeros(0,num_sec+3);
i_3=zeros(0,num_sec+3);
i_4=zeros(0,num_sec+3);
i_5=zeros(0,num_sec+3);
i_6=zeros(0,num_sec+3);
i_7=zeros(0,num_sec+3);

%%A=semi-major axis
A=[1, 5263138304; 2, 9832684544; 3, 13599840256; 4, 20726155264; 5, 40839348203; 6, 68773560320; 7, 90118820000];
%%e_e=eccentricity
e_e=[1, .2; 2, .01; 3, 0; 4, .05; 5, .14; 6, .05; 7, .26];
%%I_i=inclination
I_i=[1,7; 2, 2.1; 3, 0; 4, .06; 5, 5; 6, 1.304; 7, 6.15];
%%m_0=mean anomoly
m_0=[1, pi; 2, pi; 3, pi; 4, pi; 5, pi; 6, .1; 7, pi];
%%argument of periapse
A_p=[1, 15*pi/180; 2, 0; 3, 0; 4, 0; 5, .5*pi; 6, 0; 7, 260*pi/180];
%%Longitude of Ascending Node
LAN=[1, 70*pi/180; 2, 15*pi/180; 3, 0; 4, 135.5*pi/180; 5, 80*pi/180; 6, 52*pi/180; 7, 260*pi/180];


%%semi-major axis a
for planet=3%1:1:7
    t=0;
    
    
    % planet = input('what planet?     ');
    a=A(planet,2);
    e=e_e(planet,2);
    I=I_i(planet,2);
    %   a=input('What is the semi-major axis?     ');%%Semi-Major Axis
    %   e=input('what is the eccentricity?     ');%%Eccentricity
    %   I=input('what is the inclination?     ');%%Incliination
    
    
    
    count=0;
    %%i = tan(I+(pi/2));
    I=I.*(2*pi/360);
    interval_time=clock;
    planet
    clock_time=etime(interval_time,clockstart)
    
    
    for t=0:160000000/num_sec:160000000
        t
        syms E;
        n=sqrt((G*M)/(a.^3));
        count=count+1;
        T(count)=t;
        m=mod(n*(t)+m_0(planet,2), 2.*pi);
        E=solve(m==E-e*sin(E));
        %%v=2*atan(sqrt(1-e)*cos(E/2)/(sqrt(1+e)*sin(E/2)));
        %%v=acos((cos(E-e))/(1+e*cos(E)));
        v=mod(double(2*atan(sqrt((1+e)/(1-e))*tan(E/2))+2*pi) ,2*pi);
        
        r=(a*(1-e.^2))/(1+e.*cos(v));
        if planet==1
            v_1(count)=v;
            R_1(count)=r;
            ang_1(count)=mod(v+LAN(planet,2),2.*pi);
            
            i_mod=-I*sin(v+A_p(planet,2));
            i_1(count)=i_mod;
            H=r*sin(i_mod);
            h_1(count)=H;
        elseif planet==2
            v_2(count)=v;
            R_2(count)=r;
            ang_2(count)=mod(v+LAN(planet,2),2.*pi);
            i_mod=-I*sin(v+A_p(planet,2));
            i_2(count)=i_mod;
            H=r*sin(i_mod);
            h_2(count)=H;
        elseif planet==3
            v_3(count)=v;
            R_3(count)=r;
            ang_3(count)=mod(v+LAN(planet,2),2.*pi);
            i_mod=-I*sin(v+A_p(planet,2));
            i_3(count)=i_mod;
            H=r*sin(i_mod);
            h_3(count)=H;
%         elseif planet==4
%             v_4(count)=v;
%             R_4(count)=r;
%             ang_4(count)=mod(v+LAN(planet,2),2.*pi);
%             i_mod=-I*sin(v+A_p(planet,2));
%             i_4(count)=i_mod;
%             H=r*sin(i_mod);
%             h_4(count)=H;
%         elseif planet==5
%             v_5(count)=v;
%             R_5(count)=r;
%             ang_5(count)=mod(v+LAN(planet,2),2.*pi);
%             i_mod=-I*sin(v+A_p(planet,2));
%             i_5(count)=i_mod;
%             H=r*sin(i_mod);
%             h_5(count)=H;
%         elseif planet==6
%             v_6(count)=v;
%             R_6(count)=r;
%             ang_6(count)=mod(v+LAN(planet,2),2.*pi);
%             i_mod=-I*sin(v+A_p(planet,2));
%             i_6(count)=i_mod;
%             H=r*sin(i_mod);
%             h_6(count)=H;
%         elseif planet==7
%             v_7(count)=v;
%             R_7(count)=r;
%             ang_7(count)=mod(v+LAN(planet,2),2.*pi);
%             i_mod=-I*sin(v+A_p(planet,2));
%             i_7(count)=i_mod;
%             H=r*sin(i_mod);
%             h_7(count)=H;
        end
        
    end
end
T = transpose(T);

% consolidated_1 = [T transpose(R_1) transpose(ang_1) transpose(h_1)];
% consolidated_2 = [T transpose(R_2) transpose(ang_2) transpose(h_2)];
consolidated_3 = [T transpose(R_3) transpose(ang_3) transpose(h_3)];
%consolidated_4 = [T transpose(R_4) transpose(ang_4) transpose(h_4)];
% consolidated_5 = [T transpose(R_5) transpose(ang_5) transpose(h_5)];
% consolidated_6 = [T transpose(R_6) transpose(ang_6) transpose(h_6)];
% consolidated_7 = [T transpose(R_7) transpose(ang_7) transpose(h_7)];
% [X_1,Y_1,Z_1]=pol2cart(ang_1,R_1,h_1);
% [X_2,Y_2,Z_2]=pol2cart(ang_2,R_2,h_2);
% [X_3,Y_3,Z_3]=pol2cart(ang_3,R_3,h_3);
% [X_4,Y_4,Z_4]=pol2cart(ang_4,R_4,h_4);
% [X_5,Y_5,Z_5]=pol2cart(ang_5,R_5,h_5);
% [X_6,Y_6,Z_6]=pol2cart(ang_6,R_6,h_6);
% [X_7,Y_7,Z_7]=pol2cart(ang_7,R_7,h_7);
% 
% 
% plot3(X_1,Y_1,Z_1)
% hold all
% plot3(X_2,Y_2,Z_2)
% plot3(X_3,Y_3,Z_3)
% plot3(X_4,Y_4,Z_4)
% plot3(X_5,Y_5,Z_5)
% plot3(X_6,Y_6,Z_6)
% plot3(X_7,Y_7,Z_7)
% axis([1.5*-10^11 1.5*10^11 1.5*-10^11 1.5*10^11 1.5*-10^11 1.5*10^11])
%%plot(T,ang)
endtime = clock
time_elapsed = etime(endtime,clockstart)
