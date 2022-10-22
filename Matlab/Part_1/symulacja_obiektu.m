clc;
clear all;

%Constants
A1 = 505;
C2 = 606;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy
tau = 120;
h2_0 = 38.44;
F1 = 78;
FD = 15;

%Obliczanie Dynamiki Objętości
funkcja = @(t,h)stan_ciagly(t,h,F1,FD,A1,C2,ap1,ap2,tau) ;
[t,h] = ode45(funkcja,[0, 5000],[0;h2_0]);

%Plot h/t

figure;
plot(t,h(:,1));
hold on;
plot(t,h(:,2));
legend;