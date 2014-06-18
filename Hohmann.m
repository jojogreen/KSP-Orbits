R1 = 13599840256;%Radius of Kerbin (and spaceship
Planet4 = importdata('consolidated_4.mat');
% Planet2 = importdata(R2);
 Planet3 = importdata('consolidated_3.mat');
% Planet4 = importdata(R4);
% Planet5 = importdata(R5);
% Planet6 = importdata(R6);
% Planet7 = importdata(R7);
mu = 1.17*10^18;
R2guess = [1.969*10^10:20000:2.1762*10^10];
T = Planet4(:,1);
R = Planet4(:,2);
ang = Planet4(:,3);
h = Planet4(:,4);
angShip =Planet3(:,3);

Rad = sqrt(h.^2 + R.^2);
Lookup = [T Rad];
RSSS = [];
deltV = [];
Tinit = [];
deltT = pi*sqrt(((R1+R2guess).^3)./(8*mu));
PHASE = [];
diffang = mod((ang-angShip), 2.*pi);
for Ti = 3*10^6:20000:5*10^6
        Tf = Ti+deltT;
        Ti
        for i =1:10360
            tf = Tf(1,i);
            [row col] = find(tf+55>T & tf-55<T);
            if abs(T(row)-tf)<=55
                
                R2real = Lookup(row,2);
                phase = pi-sqrt(mu/R2real)*(tf-Ti)/R2real;
                PHASE = [PHASE; phase];
                if abs(diffang)>phase+5*pi/180 & abs(diffang)<phase-5*pi/180
                    RSSS=[RSSS; R2real];
                    Tinit = [Tinit;Ti];
                    dV = sqrt(mu/R1)*(sqrt((2*R2real)/(R1+R2real))-1);
                    deltV = [deltV; dV];
                end
            end
        end
end
coagulated = [Tinit RSSS deltV];
scatter3(Tinit, RSSS, deltV, '.');
xlabel('Time of burn');
ylabel('Radius of arrival');
zlabel('deltaV');
        
        

%end






% 
% deltT = pi*sqrt(((R1+Rad2_4).^3)./(8*mu));
% Ti = T-deltT;
% 
% deltV = sqrt(mu/R1)*(sqrt((2*Rad2_4)./(R1+Rad2_4)) -1);
