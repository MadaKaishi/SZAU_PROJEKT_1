clear; clc;

%Parametry regulatora
Nu = 450;
N = 600;
D = N;
lamb = 100;

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

t_sym = 3000; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = 120/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:T:t_sym/T) = F1;
FD(1:T:t_sym/T) = FD;


%Wartosc zadana
yzad(1:ks)=78; yzad(ks:kk)=80;


%Symulacja obiektu
for k = kp:t_sym/T
    F1in(k) = F1+dF1;
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FD(k-1) - ap1*sqrt(v1(k-1)/A1));
    v2(k) = v2(k-1) + T*(ap1*sqrt(v1(k-1)/A1) - ap2*(nthroot(v2(k-1)/C2,4)));  
end
%Obliczanie Wysokości na podstawie objętości
h2(:,1) = sqrt(v2 /C2);

%Skok jednostkowy
Yj = h2 - min(h2);
Yj = Yj/dF1;
s = Yj(1:D,:);
plot(Yj)

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

%warunki_początkowe
kp = 300/T + 2;
ks = max(19,D+7); %chwila skoku wartosci zadania
kk = t_sym/T;
v1(1:kp) = v1_0;
v2(1:kp) = v2_0;
F1in(1:1000/T) = F1;
F1in(1000/T:kk) = F1+2;
FD(1:T:t_sym/T) = FD;


%główne wykonanie programu
for k=kp:kk
    for n=1:N
    %yzad dla horyzontu predykcji
        Y_zad(n,1) = yzad(k);
    end
    %symulacja obiektu
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FD(k-1) - ap1*sqrt(v1(k-1)/A1));
    v2(k) = v2(k-1) + T*(ap1*sqrt(v1(k-1)/A1) - ap2*(sqrt(sqrt(v2(k-1)/C2)))); 
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
    wejscie_F1in(k)=F1in(k);
    wyjscie_h2(k)=h2(k); 
end


iteracja = 0:1:kk-1;
%Plot wyjście
figure;
stairs(iteracja, wyjscie_h2)
hold on;
stairs(iteracja, yzad);
hold off;
xlabel('k'); ylabel("y");
