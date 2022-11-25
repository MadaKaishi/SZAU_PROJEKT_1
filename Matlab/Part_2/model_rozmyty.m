clc;
clear all;

%% Parametry programu

draw = true;
sa = false;
il_fun = 5;
set(0,'DefaultStairLineWidth',1);

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
v2_0 = h2_0^2 * C2;
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
hr01 = hr0.*m;
vr2 = hr0.^2 * C2;
Fr0 = ap1*hr01.^0.5-FD0;

%% Symulacja

if draw
    figure
    title("Porownanie obiektów")
    xlabel('k')
    ylabel('h_2')
    hold on
end

kp = tau/T + 2;
kk = t_sym/T;


for P = 36:21:120

    h1 = h1_0 * ones(il_fun+3,kk);
    h2 = h2_0 * ones(il_fun+3,kk);
    v1 = h1_0 * A1 * ones(il_fun+3,kk);
    v2 = h2_0^2 * C2 * ones(il_fun+3,kk);

    F1in(1:kk) = F1;
    FDc(1:kk) = FD;
    
    for k = kp:kk
        
        if k/T > 180
            F1in(k) = P;
        end

        % Model nieliniowy
        v1(il_fun+2,k) = v1(il_fun+2,k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*(sqrt(h1(il_fun+2,k-1))));
        v2(il_fun+2,k) = v2(il_fun+2,k-1) + T*(ap1*sqrt(h1(il_fun+2,k-1)) - ap2*(sqrt(h2(il_fun+2,k-1)))); 
        h1(il_fun+2,k) = v1(il_fun+2,k)/A1;
        h2(il_fun+2,k) = sqrt(v2(il_fun+2,k)/C2);

        % Model liniowy
        v1(il_fun+3,k) = v1(il_fun+3,k-1) + T*(F1in(k-1-(tau/T)) - F10 + FDc(k-1) - FD0 - (ap1/(2*(sqrt(h1_0))))*(h1(il_fun+3,k-1)-h1_0));
        v2(il_fun+3,k) = v2(il_fun+3,k-1) + T*((ap1/(2*(sqrt(h1_0))))*(h1(il_fun+3,k-1)-h1_0) - (ap2/(2*(sqrt(h2_0))))*(h2(il_fun+3,k-1)-h2_0)); 
        h2(il_fun+3,k) = h2_0 + 1/(2*sqrt(C2*v2_0))*(v2(il_fun+3,k) - v2_0);
        h1(il_fun+3,k) = v1(il_fun+3,k)/A1;

        for i = 1:il_fun
            %Rownania modelu
       
            v1(i,k) = v1(il_fun+1,k-1) + T*(F1in(k-1-(tau/T)) - Fr0(i) + FDc(k-1) - FD0 - (ap1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)));
            v2(i,k) = v2(il_fun+1,k-1) + T*((ap1/(2*(sqrt(hr01(i)))))*(h1(il_fun+1,k-1)-hr01(i)) - (ap2/(2*(sqrt(hr0(i)))))*(h2(il_fun+1,k-1)-hr0(i))); 
            h2(i,k) = hr0(i) + (v2(i,k) - vr2(i))*1/(2*sqrt(C2*vr2(i)));
            h1(i,k) = v1(i,k)/A1;

            %Liczenie funkcji przynaleznosci
            if i == 1
                w(i) = trapmf(h2(i,k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
            elseif i == il_fun
                w(i) = trapmf(h2(i,k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
            else
                w(i) = trapmf(h2(i,k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
            end
        end
        %Wyliczanie wyjscia modelu

        h2(il_fun+1,k) = w * h2(1:il_fun, k)/sum(w);
        h1(il_fun+1,k) = w * h1(1:il_fun, k)/sum(w);
        v2(il_fun+1,k) = w * v2(1:il_fun, k)/sum(w);
        v1(il_fun+1,k) = w * v1(1:il_fun, k)/sum(w);

    end


    stairs((1:k),h2(il_fun+2,:),"--r")
    hold on
    stairs((1:k),h2(il_fun+3,:),"--g")
    stairs((1:k),h2(il_fun+1,:),"b")
    legend("Model nieliniowy", "Model zlinearyzowany", "Model rozmyty","Location","northoutside","Orientation","horizontal")

    clear v1 v2 h1 h2
end
if sa
    print(sprintf("Porown_rozm_mod_il_%i",il_fun),'-dpng','-r400');
end