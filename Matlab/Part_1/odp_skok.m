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
v2_0 = h2_0^2 * C2;
v1_0 = h1_0 * A1;

t_sym = 1500; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:T:t_sym/T) = F1;
FDc(1:T:t_sym/T) = FD;


%Symulacja obiektu
for k = kp:t_sym/T
    F1in(k) = F1+dF1;
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*sqrt(v1(k-1)/A1));
    v2(k) = v2(k-1) + T*(ap1*sqrt(v1(k-1)/A1) - ap2*(nthroot(v2(k-1)/C2,4)));  
end
%Obliczanie Wysokości na podstawie objętości
h2(:,1) = sqrt(v2 /C2);

%Skok jednostkowy
Yj = h2 - min(h2);
Yj = Yj/dF1;
s = Yj(1:D,:);
plot(Yj)

save("Odp_skok\odp_skok.mat", "s");