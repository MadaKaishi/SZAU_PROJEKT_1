clear all

%% Wymierzanie skoku jednostkowego

D = 2000;

%skok wartosci F1
dF1 = 1;

%Parametry obiektu
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy obiektu
tau = 120;
h2_0 = 38.44;
h1_0 = 16.34;
F1 = 78;
FD = 15;

F10 = 78;
FD0 = 15;

%Objetośc punktu pracy
v2_0 = h2_0^2 * C2;
v1_0 = h1_0 * A1;

t_sym = 2200; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
kk = t_sym/T;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:T:t_sym/T) = F1;
FDc(1:T:t_sym/T) = FD;


%Symulacja obiektu liniowego
for k = kp:t_sym/T
    F1in(k) = F1+dF1;
        v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) - F10 + FDc(k-1) - FD0 - (ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0));
        v2(k) = v2(k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(k-1)-h2_0)); 
        h2(k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(k) - v2_0);
        h1(k) = v1(k)/A1;
end


%Skok jednostkowy
Yj = h2 - h2(kp);
Yj = Yj/dF1;
s = Yj(:,kp:D+1);
plot(s)
title("Odpowiedz skokowa")
ylabel("h_2")
xlabel("k")
set(get(gca,'ylabel'),'rotation',0)
xlim([0 1800])
exportgraphics(gca,'odp_skok_lin.pdf')

% save("Odp_skok\odp_skok.mat", "s");