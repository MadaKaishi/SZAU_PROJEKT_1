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
kp = tau/T + 2;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
F1in(1:T:t_sym/T) = F1;
FDc(1:T:t_sym/T) = FD;

for k = kp:t_sym/T
    if k/T > 180
        F1in(k) = 78;
    end
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) - F1 + FDc(k-1) - FD - (ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0));
    v2(k) = v2(k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(k-1)-h2_0)); 
    h2(k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(k) - v2_0);
    h1(k) = v1(k)/A1;

end

%Obliczanie Wysokości na podstawie objętości
h(:,1) = v1 / A1;
h(:,2) = sqrt(v2 /C2);

figure;
set(0,'defaultLineLineWidth',1);
set(0,'DefaultStairLineWidth',1);
plot((1:k),h(:,1));
hold on;
plot((1:k),h(:,2));
legend("h_1","h_2");
xlabel("k"); ylabel("Wysokość");
%print('h1_h2_niel.png','-dpng','-r400')