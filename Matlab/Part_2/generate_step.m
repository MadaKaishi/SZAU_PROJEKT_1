function [s] = generate_step(Y0, draw)

%% Wymierzanie skoku jednostkowego


D = 2000;
D = D + 120;

% Y0 = 120;
% draw = true;

%% Parametry obiektu
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%% Nowy punkt pracy obiektu
tau = 120;
h2_0 = Y0;
h1_0 = h2_0 * (ap2/ap1)^2;
v1_0 = h1_0 * A1;
v2_0 = h2_0^2 * C2;
FD0 = 15;
Fr0 = ap1*h1_0^0.5 - FD0;

%% Punkt pracy z miejsca skoku

t_sym = 4000; %czas symulacji
T = 1; %krok

%% Parametry stanu ustalonego
s_h2_0 = 38.44;
s_h1_0 = 16.34;
s_F10 = 78;
s_FD0 = 15;
s_v2_0 = s_h2_0^2 * C2;
s_v1_0 = s_h1_0 * A1;

%% Warunki_początkowe
kp = 120/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
k_skok = 2000;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:T:t_sym/T) = Fr0;
FDc(1:T:t_sym/T) = FD0;




%% Symulacja obiektu w punkcie pracy
for k = kp:t_sym/T
    if k >= k_skok
        F1in(k) = Fr0+1;
    end
        v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) - Fr0 + FDc(k-1) - FD0 - (ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0));
        v2(k) = v2(k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(k-1)-h2_0)); 
        h2(k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(k) - v2_0);
        h1(k) = v1(k)/A1;
end

%% Skok jednostkowy

Yj = h2 - h2(kp);
Yback = Yj;
Yj = Yj/(1);

s = Yj(:,k_skok:kk);
s_untimmed = s;
s = s - ones(1,length(s))*s(1);

if draw
    plot(s)
    title("Odpowiedz skokowa")
    ylabel("h_2")
    xlabel("k")
    set(get(gca,'ylabel'),'rotation',0)
end

end