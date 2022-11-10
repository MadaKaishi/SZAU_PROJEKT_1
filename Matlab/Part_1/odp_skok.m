clear all

%% Wymierzanie skoku jednostkowego

D = 1500;

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

%Objetośc punktu pracy
% v2_0 = h2_0^2 * C2;
% v1_0 = h1_0 * A1;

t_sym = 1600; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
F1in(1:T:t_sym/T) = F1;
FDc(1:T:t_sym/T) = FD;


%Symulacja obiektu
for k = kp:t_sym/T
    F1in(k) = F1+dF1;
        h1(k) = h1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1 * (sqrt(h1_0) + (h1(k-1) - h1_0)/(2*sqrt(h1_0))))/A1;
        h2(k) = h2(k-1) + T*((ap1 * sqrt(h1_0) - ap2 * sqrt(h2_0))/(2*C2*h2_0)   +    ap1*(h1(k-1) - h1_0)/(4*C2*(sqrt(h1_0)))    -   ap2*(h2(k-1) - h2_0)/(4*C2*(nthroot((h2_0 * h2_0),3))));
end
%Obliczanie Wysokości na podstawie objętości
% h2(:,1) = sqrt(v2 /C2);

%Skok jednostkowy
Yj = h2 - h2(1);
Yj = Yj/dF1;
s = Yj(:,1:D);
plot(Yj)
title("Odpowiedz skokowa")
ylabel("h_2")
xlabel("k")
set(get(gca,'ylabel'),'rotation',0)
exportgraphics(gca,'odp_skok_lin.pdf')

% save("Odp_skok\odp_skok.mat", "s");