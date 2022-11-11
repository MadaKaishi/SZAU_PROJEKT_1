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

%Wybranie punktu pracy
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
    title('Przebiegi wyjœcia dla skoku wartoœci sterowania w chwili 180.')
    xlabel('k')
    ylabel('h_2')
    hold on
end

kp = tau/T + 2;
kp = t_sym/t;


for P = 36:21:120
    h1 = h1_0 * ones(il_fun,kp);
    h2 = h2_0 * ones(il_fun,kp);
    F1in(1:kk) = F1;
    FDc(1:kk) = FD;
    for k = kp:kk
        if k/T > 180
            F1in(k) = P;
        end
        for i = 1:il_fun
            h1(i,k) = h1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1 * (sqrt((hr0(r)*m)) + (h1(k-1) - h1_0)/(2*sqrt(hr0(r)*m))))/A1;
            h2(i,k) = h2(k-1) + T*((ap1 * sqrt((hr0(r)*m)) - ap2 * sqrt((hr0(r))))/(2*C2*h2_0)   +    ap1*(h1(k-1) - (hr0(r)*m))/(4*C2*(sqrt((hr0(r)*m))))    -  ap2*(h2(k-1) - hr0(r))/(4*C2*(nthroot((hr0(r) * hr0(r)),3))));
            if i == 1
                w(r) = trapmf(h,[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == il_fun
                plot(trapmf(h,[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]));
            else
                plot(trapmf(h,[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]));
            end
        end
    end

end













%warunki_początkowe
kp = tau/T + 2;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
F1in(1:T:t_sym/T) = F1;
FDc(1:T:t_sym/T) = FD;
figure;

for P = 36:21:120
    h2(1:kp) = h2_0;
    h1(1:kp) = h1_0;
    for k = kp:t_sym/T
        if k/T > 180
            F1in(k) = P;
        end
%         v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) - F10 + FDc(k-1) - FD0 - (ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0));
%         v2(k) = v2(k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(k-1)-h2_0)); 
%         h2(k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(k) - v2_0);
%         h1(k) = v1(k)/A1;

        % Version 2.0
          h1(k) = h1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1 * (sqrt(h1_0) + (h1(k-1) - h1_0)/(2*sqrt(h1_0))))/A1;
          h2(k) = h2(k-1) + T*((ap1 * sqrt(h1_0) - ap2 * sqrt(h2_0))/(2*C2*h2_0)   +    ap1*(h1(k-1) - h1_0)/(4*C2*(sqrt(h1_0)))    -   ap2*(h2(k-1) - h2_0)/(4*C2*(nthroot((h2_0 * h2_0),3))));
    end
    plot((1:k),h2);
    hold on
    clear h2
    clear h1
end

xlabel("k"); ylabel("h_2");
%print('h1_h2_niel.png','-dpng','-r400')