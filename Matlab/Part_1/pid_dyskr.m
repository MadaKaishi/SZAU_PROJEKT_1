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

t_sym = 5000; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
F1in(1:kp) = F1;
FD(1:T:t_sym/T) = FD;
e(1:kp) = 0;

yzad(1:kp) = 38.44;
yzad(kp:kk) = 50;

%Nastawy PID
Kp = 1.93;
Ti = 250;
Td = 60;

r0 = Kp*(1 + T/(2*Ti) + Td/T);
r1 = Kp*(T/(2*Ti) - (2*Td)/T -1);
r2 = Kp*Td/T;

for k = kp:t_sym/T

    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FD(k-1) - ap1*sqrt(h1(k-1)));
    v2(k) = v2(k-1) + T*(ap1*sqrt(h1(k-1)) - ap2*(sqrt(h2(k-1))));
    h1(k) = v1(k)/A1;
    h2(k) = sqrt(v2(k)/C2);

    %PID
    e(k) = yzad(k) - h2(k);
    F1in(k) = F1in(k-1) + r0*e(k) + r1*e(k-1) + r2*e(k-2);

end

iteracja = 0:1:kk-1;
%Plot wyjście
figure;
stairs(iteracja, h2)
hold on;
stairs(iteracja, yzad,"--");
hold off;
xlabel('k'); ylabel("y");