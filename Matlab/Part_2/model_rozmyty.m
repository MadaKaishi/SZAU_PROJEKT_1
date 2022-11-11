clc;
clear all;

%% Parametry programu

draw = true;
sa = true;
il_fun = 2;

%% Definicja parametrów

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

F10 = 78;
FD0 = 15;

t_sym = 5000; %czas symulacji
T = 1; %krok


%% Zmienne modelu rozmytego
h_min = 0;
h_max = 90;
h = (h_min:1:h_max)';

nach = 3; %nachylenie funkcji 

d = (h_max-h_min)/il_fun; %szerokości funkcji przynależnośći
c = h_min+d:d:h_max-d; %punkty przegięcia

%Wybranie punktu linearyzacji
hr0 = ones(1,il_fun);
hr0(1) = d/2;
hr0(il_fun) = min((h_max+c(il_fun-1))/2+1, h_max);
    if il_fun > 2
        hr0(2:il_fun-1) = (c(2:il_fun-1)+c(1:il_fun-2))./2;
    end

m = (ap2/ap1)^2;


%% Warunki początkowe

if draw
    figure
    title('Przebiegi dla skoku sterowania w chwili 180.')
    xlabel('k')
    ylabel('h_2')
    hold on
end

kp = tau/T + 2;
kk = t_sym/T;


for P = 78
    h1 = h1_0 * ones(il_fun+1,kp);
    h2 = h2_0 * ones(il_fun+1,kp);
    F1in(1:kk) = F1;
    FDc(1:kk) = FD;
    for k = kp:kk
        if k/T > 180
            F1in(k) = P;
        end
        for i = 1:il_fun
            h1(i,k) = h1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1 * (sqrt((hr0(i)*m)) + (h1(k-1) - h1_0)/(2*sqrt(hr0(i)*m))))/A1;
            h2(i,k) = h2(k-1) + T*((ap1 * sqrt((hr0(i)*m)) - ap2 * sqrt((hr0(i))))/(2*C2*h2_0)   +    ap1*(h1(k-1) - (hr0(i)*m))/(4*C2*(sqrt((hr0(i)*m))))    -  ap2*(h2(k-1) - hr0(i))/(4*C2*(nthroot((hr0(i) * hr0(i)),3))));
            if i == 1
                w(i) = trapmf(1,[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == il_fun
                w(i) = trapmf(1,[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
            else
                w(i) = trapmf(1,[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
            end
        end
        h2(il_fun+1,k) = w*h2(1:il_fun, k)/sum(w);
        h1(il_fun+1,k) = w*h1(1:il_fun, k)/sum(w);
    end
end

plot(h2(3,:))
plot(h1(3,:))