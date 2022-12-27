clear; clc;

file = load("Odp_skok\odp_skok.mat");
s = file.s;

%Parametry regulatora
Nu = 1400;
N = 1400;
D = 1400;
lamb = 5;

Umax = 300;
Umin = 0;

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



%% Macierze do DMC

%Wyznaczanie macierzy M
M=zeros(N,Nu);
for i=1:N
   for j=1:Nu
      if (i>=j)
         M(i,j)=s(i-j+1);
      end
   end
end

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



K = ((M'*M + lamb * eye(Nu))^(-1))* M';
DUp = zeros(1,D-1);
Y = zeros(N,1);
Ku = K(1,:)*MP;
Ke = sum(K(1,:));

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
FDc(1:kp) = FD;
FDc(kp:5000) = 30;
FDc(5000:10000) = 15;
FDc(10000:15000) = 7.5;


yzad(1:kk)=38.44; 


error = 0;
err = 0;
del_u = 0;

%główne wykonanie programu
for k=kp:kk

    %symulacja obiektu
    v1(k) = v1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1*sqrt(h1(k-1)));
    v2(k) = v2(k-1) + T*(ap1*sqrt(h1(k-1)) - ap2*(sqrt(h2(k-1))));
    h1(k) = v1(k)/A1;
    h2(k) = sqrt(v2(k)/C2);
    
    %stała trajektoria referencyjna
    for i=D-1:-1:2
        DUp(i) = DUp(i-1);
    end
    
    DUp(1) = del_u;

    err = yzad(k) - h2(k);
    del_u = Ke*err-Ku*DUp';
   
   
    F1in(k)=F1in(k-1)+del_u;  
    error = error + norm((yzad(k) - h2(k)))^2;

    if F1in(k) > Umax
        F1in(k) = Umax;
    elseif F1in(k) < Umin
        F1in(k) = Umin;
    end

    del_u =  F1in(k) -  F1in(k-1);

end
f = figure;
subplot(3,1,1)
stairs(1:kk,h2)
xlabel("k")
ylabel("h_2")
title(sprintf("h_2 error=%d",error))
ylim([25, 50])


subplot(3,1,2)
stairs(1:kk,FDc)
xlabel("k")
ylabel("Fd")
title("Fd")
ylim([0, 32])


subplot(3,1,3)
stairs(1:kk,F1in)
xlabel("k")
ylabel("F_1_i_n")
title("F_1_i_n")
ylim([60, 90])

display(error)