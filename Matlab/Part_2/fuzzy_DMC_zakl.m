clear; clc;
%% Parametry programu
draw = true;
sa = false;
draw_f_przyn = true;

set(0,'DefaultStairLineWidth',1);
DUmax = 1;
Umax = 120;
Umin = 0;

%% Parametry regulatora
Nu = 4;
N = 1200;
D = 1500;
lamb = 100;

%liczba regulatorów
il_fun = 5;
lamb = lamb*ones(1,il_fun);


%% Parametry modelu i symulacji

%Constants
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

%Podstawowe punkty linearyzacji
F10 = 78;
FD0 = 15;

%Objetośc punktu pracy
v2_0 = h2_0^2 * C2;
v1_0 = h1_0 * A1;

t_sym = 15000; %czas symulacji
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


ku = zeros(il_fun,D-1);
ke = zeros(1,il_fun);

%% Pokazanie funkcji przynależności

if draw_f_przyn
    figure
    hold on
    %Plotter funkcji przynaleznosci
    for i = 1:il_fun
        if i == 1
            plot(trapmf(h_min:1:h_max,[0 0 c(1)-nach/2 c(1)+ nach/2]));
        elseif i == il_fun
            plot(trapmf(h_min:1:h_max,[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]));
        else
            plot(trapmf(h_min:1:h_max,[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]));
        end
    end
    xlim([0 90])
    xlabel("h_2"); ylabel("Funkcja przynależności");
    title(sprintf("Funkcja przynaleznosci dla %i zbiorów rozmytych",il_fun))
    if sa
        print(sprintf('funkcja_przynelznosci_%i.png',il_fun),'-dpng','-r400')
    end
end


%% Liczenie poszczególnych regulatorów

for r = 1:il_fun
    s = generate_step(hr0(r),false);
        
    M=zeros(N,Nu);
        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j)=s(i-j+1);
              end
           end
        end
        

    MP=zeros(N,D-1);
        for i=1:N
           for j=1:D-1
              if i+j<=D
                 MP(i,j)=s(i+j)-s(j);
              else
                 MP(i,j)=s(D)-s(j);
              end      
           end
        end
    
    K = ((M'*M + lamb(r) * eye(Nu))^(-1))* M';
    ku(r,:) = K(1,:)*MP;
    ke(r) = sum(K(1,:));
end

%% Symulacja obiektu

%warunki_początkowe
kp = D/T + 2;
ks = max(19,D+100); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h2(1:kp) = h2_0;
h1(1:kp) = h1_0;
F1in(1:T:kp) = F1;
FD = 15;
FDc(kp:5000) = 30;
FDc(5000:10000) = 15;
FDc(10000:15000) = 7.5;

%Skok wartosci zadanej:
yzad(1:kk)=38.44; 



error = 0;
w = zeros(1,il_fun);
Du = zeros(il_fun,1);
DUp = zeros(1,D-1);
Y = zeros(N,1);
err_cur = 0;
err_sum = 0;

%główne wykonanie programu
for k=kp:kk
    for n=1:N


    end
    %symulacja obiektu
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*sqrt(h1(k-1)));
    v2(k) = v2(k-1) + T*(ap1*sqrt(h1(k-1)) - ap2*(sqrt(h2(k-1))));
    h1(k) = v1(k)/A1;
    h2(k) = sqrt(v2(k)/C2);
    
    %Liczenie błędu
    err_cur = yzad(k) - h2(k);
    err_sum = err_sum + norm((yzad(k) - h2(k)))^2;


    %Liczenie wartości przyrostu sterowania
    for i = 1:il_fun

        Du(i) = ke(i)*err_cur-ku(i,:)*DUp';

        if i == 1
            w(i) = trapmf(h2(k),[0 0 c(1)-nach/2 c(1)+ nach/2]);
        elseif i == il_fun
            w(i) = trapmf(h2(k),[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]);
        else
            w(i) = trapmf(h2(k),[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]);
        end
    end
    
    %Sprawdzenie liczenia wag
    w_over_time(:,k) = w;

    %Ogranieczenia przyrostu sterowania
    DUfin = w * Du / sum(w);
    
    if DUfin > DUmax
        DUfin = DUmax;
    elseif DUfin < -DUmax
        DUfin = -DUmax;
    end

    for i = D-1:-1:2
      DUp(i) = DUp(i-1);
    end
    DUp(1) = DUfin;

    F1in(k) = F1in(k-1) + DUp(1);

    %Ograniczenia sterowania
    if F1in(k) > Umax
        F1in(k) = Umax;
    elseif F1in(k) < Umin
        F1in(k) = Umin;
    end

end

if draw
%Plot wyjście
f = figure;
subplot(3,1,1)
stairs(1:kk,h2)
xlabel("k")
ylabel("h_2")
title(sprintf("h_2 error=%d",err_sum))
ylim([25, 50])

subplot(3,1,2)
stairs(1:kk,FDc)
xlabel("k")
ylabel("Fd")
title("Fd")
ylim([0, 32])
% exportgraphics(f,'odp_na_zmiane_zakl.pdf')

subplot(3,1,3)
stairs(1:kk,F1in)
xlabel("k")
ylabel("F_1_i_n")
title("F_1_i_n")
ylim([60, 90])

% exportgraphics(f,'odp_na_zmiane_zakl.pdf')

end
display(err_sum)

