clc;
clear all;

%Constants
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy
tau = 120;
h2_0 = 38.44;
h1_0 = 16.34;
F1 = 78;
FD = 15;

%Objetośc punktu pracy
v2_0 = h2_0^2 * C2;
v1_0 = h1_0 * A1;

for F1 = 36:21:120
%Obliczanie Dynamiki Objętości
funkcja = @(t,v)stan_ciagly(t,v,F1,FD,A1,C2,ap1,ap2,tau) ;
[t,v] = ode45(funkcja,[0, 2000],[v1_0;v2_0]);

%Obliczanie Wysokości na podstawie objętości
h(:,2) = sqrt(v(:,2) /C2);

%Plot h/t
set(0,'defaultLineLineWidth',1);
set(0,'DefaultStairLineWidth',1);
plot(t,h(:,2));
hold on
name = "F_1 = "+ F1;
clear h;
end

xlabel("Czas"); ylabel("h_2");

lgd = cell(5,1) ;
for i=1:5
    lgd{i} = strcat('F_1= ',num2str(15+21*i)) ;
end
legend(lgd,Location="northoutside",Orientation="horizontal")
xlim([0 1500])
% print('h2_niel.png','-dpng','-r400')