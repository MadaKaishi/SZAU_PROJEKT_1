function [s] = generate_step(Y0, draw)

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
h2_0 = Y0;
h1_0 = h2_0 * (ap2/ap1)^2;
v1_0 = h1_0 * A1;
v2_0 = h2_0^2 * C2;
FD0 = 15;
F1 = 78;
Fr0 = ap1*h1_0^0.5 - FD0;


t_sym = 2000; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:T:t_sym/T) = Fr0;
FDc(1:T:t_sym/T) = FD0;


%Symulacja obiektu w punkcie pracy
for k = kp:t_sym/T
    F1in(k) = Fr0+dF1;
        v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) - Fr0 + FDc(k-1) - FD0 - (ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0));
        v2(k) = v2(k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(k-1)-h2_0)); 
        h2(k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(k) - v2_0);
        h1(k) = v1(k)/A1;
end
%Obliczanie Wysokości na podstawie objętości
% h2(:,1) = sqrt(v2 /C2);

%Skok jednostkowy
Yj = h2 - h2(1);
Yj = Yj/dF1;
s = Yj(:,1:D);

if draw
    plot(Yj)
    title("Odpowiedz skokowa")
    ylabel("h_2")
    xlabel("k")
    set(get(gca,'ylabel'),'rotation',0)
end

end