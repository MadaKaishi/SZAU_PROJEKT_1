clc;
clear all;

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

%Objetośc punktu pracy
v2_0 = h2_0^2 * C2;
v1_0 = h1_0 * A1;

for F1 = 36:21:120
%Obliczanie Dynamiki Objętości

%Funkcja zlinearyzowana
funkcja = @(t_lin,h_lin)linearyzacja(t_lin,h_lin,F1,FD,A1,C2,ap1,ap2,tau,h1_0,h2_0);
[t_lin,h_lin] = ode45(funkcja,[0, 2000],[h1_0;h2_0]);

%Funkcja nieliniowa
funkcja2 = @(t_niel,v_niel)stan_ciagly(t_niel,v_niel,F1,FD,A1,C2,ap1,ap2,tau) ;
[t_niel,v_niel] = ode45(funkcja2,[0, 2000],[v1_0;v2_0]);

%Obliczanie Wysokości na podstawie objętości
h_niel(:,2) = sqrt(v_niel(:,2) /C2);

%Plot h/t

plot(t_niel,h_niel(:,2),"Color","r");
hold on
plot(t_lin,h_lin(:,2),"LineStyle","--","Color","b");
name = "F_1 = "+ F1;
clear h_lin;
clear h_niel;
end
xlabel("Czas"); ylabel("h_2");
title("Porownanie obiektów")
legend("Obiekt nieliniowy","Obiekt zlinearyzowany","Location","northoutside")
% print('porow_2_obiekt.png','-dpng','-r400')