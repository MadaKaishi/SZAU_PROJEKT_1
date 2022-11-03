clear; clc;

%Parametry regulatora
Nu = 2;
N = 8;
D = 1500;
lamb = 100;

%% Wymierzanie skoku jednostkowego

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

t_sym = 15000; %czas symulacji
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
% plot(Yj)
%print('odp_skok.png','-dpng','-r400')

%% Macierze do DMC

%Macierz MP
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

%Wyznaczanie macierzy M
M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end
   end
end

K = ((M'*M + lamb * eye(Nu))^(-1))* M';
DUp = zeros(D-1, 1);
Y = zeros(N,1);

%% Obiekt star

%warunki_początkowe
kp = D/T + 2;
ks = max(19,D+100); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
h2(1:kp) = h2_0;
h1(1:kp) = h1_0;
F1in(1:1000/T) = F1;
F1in(1000/T:kk) = F1;
FD = 15;
FDc(1:T:t_sym/T) = FD;

%Skok wartosci zadanej:
yzad(1:ks)=38.44; 
yzad(ks:5000)=50;
yzad(5000:10000)=30;
yzad(10000:15000)=40;

%główne wykonanie programu
for k=kp:kk
    for n=1:N
    %yzad dla horyzontu predykcji
        Y_zad(n,1) = yzad(k);
    end
    %symulacja obiektu
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*sqrt(h1(k-1)));
    v2(k) = v2(k-1) + T*(ap1*sqrt(h1(k-1)) - ap2*(sqrt(h2(k-1))));
    h1(k) = v1(k)/A1;
    h2(k) = sqrt(v2(k)/C2);
    %stała trajektoria referencyjna
    for n=1:N
        Y(n) = h2(k);
    end
    %DMC
    for n = 1:D-1
        DUp(n) = F1in(k-n) - F1in(k-n-1);
    end
    Yo = MP*DUp+Y;
    DU = K*(Y_zad - Yo);
    F1in(k)=F1in(k-1)+DU(1);   
end


iteracja = 0:1:kk-1;
%Plot wyjście
figure;
stairs(iteracja, h2)
hold on;
stairs(iteracja, yzad,"--");
hold off;
xlabel('k'); ylabel("h");
legend("h_2","h_2_z_a_d")
print('DMC_reg.png','-dpng','-r400')

%Plot sterowanie
figure;
stairs(iteracja, F1in)
legend("F_1_i_n")
xlabel('k'); ylabel("F_1_i_n");
print('DMC_ster.png','-dpng','-r400')