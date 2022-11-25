clear; clc;
%% Parametry programu
draw = true;
sa = false;
draw_f_przyn = false;
set(0,'DefaultStairLineWidth',1);

%% Zmienne modelu rozmytego

%liczba regulatorów
il_fun = 5;

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


%% Parametry regulatora

N = 1200;
D = 1500;
lamb = 25;
Nu = 5;


lamb = lamb*ones(1,il_fun);

MP = zeros(N,D-1,il_fun);
M = zeros(N,Nu,il_fun);

for i = 1:il_fun
    s = generate_step(hr0(i),false);
        
    for l=1:N
       for j=1:Nu
          if (l>=j)
             M(l,j,i)=s(l-j+1);
          end
       end
    end


    for l=1:N
       for j=1:D-1
          if l+j<=D
             MP(l,j,i)=s(l+j)-s(j);
          else
             MP(l,j,i)=s(D)-s(j);
          end    
       end
    end

end


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

t_sym = 20000; %czas symulacji
T = 1; %krok

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
FDc(1:T:t_sym/T) = FD;

%Skok wartosci zadanej:
yzad(1:ks)=38.44; 
yzad(ks:5000)=40;
yzad(5000:10000)=80;
yzad(10000:15000)=20;
yzad(15000:20000)=40;


error = 0;
w = zeros(1,il_fun);
Du = zeros(1,Nu)';
DUp = zeros(1,D-1)';
err_cur = 0;
err_sum = 0;
Yz = zeros(1,N)';
yk = zeros(1,N)';

A = [tril(ones(Nu));tril(ones(Nu))*-1]; %?
B = zeros(2*Nu,1); %?
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

    Mr = zeros(N,Nu);
    MPr = zeros(N,D-1);
    for i = 1:il_fun
        Mr = Mr + w(i)*M(:,:,i)/sum(w);
        MPr = MPr + w(i)*MP(:,:,i)/sum(w);
    end
    lambr = w*lamb'/sum(w);

    B(1:Nu)=(120-F1in(k-1)); %Górne ograniczenia
    B(Nu+1:end) = (F1in(k-1)-30); %Dolne ograniczenia
    Yz(1:end)=yzad(k);
    yk(1:end)=h2(k);
    Du = fmincon(@(Du)(Yz-yk-MPr*DUp-Mr*Du)'*(Yz-yk-MPr*DUp-Mr*Du)+lambr*Du'*Du,Du,A,B);
    holder = Du(1);
    for i = D-1:-1:2
        DUp(i) = DUp(i-1);
    end
    DUp(1) = holder;
    F1in(k) = F1in(k-1) + DUp(1);

    plot(h2, 'b')
    hold on
    plot(F1in,'g')
    plot(yzad,"--r")
    drawnow;

end

display(err_sum)

