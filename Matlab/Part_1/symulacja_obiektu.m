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

ode45(@dy,[0,10],[0,0])