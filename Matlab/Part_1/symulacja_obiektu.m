clc;
clear all;

%Constants
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy
F1 = 78;
FD = 15;

%Obliczanie Dynamiki Objętości
funkcja = @(t,y) [F1 + FD - ap1 * sqrt(y(1)/A1); ap1 * sqrt(y(1)/A1)- ap2*(sqrt(sqrt(y(2)/C2)))];
[t,v] = ode45(funkcja,[0, 5000],[0 0]);

%Obliczanie wysokości
h(:,1) = v(:,1)/A1;
h(:,2) = sqrt(v(:,2)/C2);

%Plot h/t

plot(t,h);