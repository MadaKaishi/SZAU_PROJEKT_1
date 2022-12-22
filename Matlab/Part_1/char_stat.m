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


set(0,'defaultLineLineWidth',1);
set(0,'DefaultStairLineWidth',1);

i = 1;
for value = 20:0.1:120
    %Warunki początkowe
    v1(1:kp) = v1_0;
    v2(1:kp) = v2_0;
    h1(1:kp) = h1_0;
    h2(1:kp) = h2_0;
    F1in(1:T:t_sym/T) = value;
    FD(1:T:t_sym/T) = FD;
    
    for k = kp:t_sym/T
        v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FD(k-1) - ap1*sqrt(h1(k-1)));
        v2(k) = v2(k-1) + T*(ap1*sqrt(h1(k-1)) - ap2*(sqrt(h2(k-1))));
        h1(k) = v1(k)/A1;
        h2(k) = sqrt(v2(k)/C2);
    end


    h2_konc(i) = h2(end);
    i = i+1;

    clear h2
    clear h1
    clear v1
    clear v2

end

plot(20:0.1:120,h2_konc)
xlabel("u")
ylabel("h_2")
title("Charakterystyka statyczna obiektu nieliniowego")