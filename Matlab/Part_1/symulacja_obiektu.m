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
h1_0 = 0;

v2_0 = h2_0^2 * C2;
v1_0 = h1_0^2 * A1;

F1 = 78;
FD = 15;

%Obliczanie Dynamiki Objętości
funkcja = @(t,v)stan_ciagly(t,v,F1,FD,A1,C2,ap1,ap2,tau) ;
[t,v] = ode45(funkcja,[0, 10000],[v1_0;v2_0]);

%Obliczanie Wysokości na podstawie objętości
h(:,1) = v(:,1) / A1;
h(:,2) = sqrt(v(:,2) /C2);

%Plot h/t
figure;
plot(t,h(:,1));
hold on;
plot(t,h(:,2));
legend("h_1","h_2");
xlabel("Czas"); ylabel("Wysokość");