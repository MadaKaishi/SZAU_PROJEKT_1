clear all

SLv = load("Pliki_mat_SL_zakl/SL_zakl.mat");

generate = false;


if generate
clear all


DMC_zakl
save("Do_porownania\DMC_zakl.mat", "h2","F1in","error");
clear all

FDMC_zakl
save("Do_porownania\FDMC_zakl.mat", "h2","F1in","err_sum");

clear all
close all

end

ks = 1500;
kp = 1402;
kk = 5000;

FD = 15;
yzad(1:5000)=38.44; 
FDc(1:ks) = FD;
FDc(ks:3250) = FD+7.5;
FDc(3250:5000) = FD-7.5;

draws = true;

porown = false;


if porown
DMCv = load("Do_porownania/DMC_zakl.mat");
FDMCv = load("Do_porownania/FDMC_zakl.mat");
SLv = load("Pliki_mat_SL_zakl/SL_zakl.mat");



f = figure;
subplot(3,1,1)
hold on;
stairs(1:kk,DMCv.h2)
stairs(1:kk, FDMCv.h2)
stairs(1:kk, SLv.h2)
stairs(yzad,"--")
hold off
xlabel("k")
ylabel("h_2")
legend("DMC","FDMC","SL",Location="west",Orientation="vertical")
title("Regulaotr DMC, error = " + DMCv.error + newline +"Regulator FDMC, error = " + FDMCv.err_sum + newline +"Regulator SL, error = " + SLv.err_sum)
ylim([25, 50])

subplot(3,1,2)
stairs(FDc)
xlabel("k")
ylabel("FD")
title("FD")
ylim([0, 37])

subplot(3,1,3)
hold on
stairs(1:kk,DMCv.F1in)
stairs(1:kk, FDMCv.F1in)
stairs(1:kk,SLv.F1in)
hold off
xlabel("k")
ylabel("F1in")
title("F1in")
legend("DMC","FDMC","SL",Location="south",Orientation="horizontal")
ylim([50, 90])
exportgraphics(f,'SL_zakl_por.pdf')

else
SLv = load("Pliki_mat_SL_zakl/SL_zakl.mat");
%Plot wyjście
f = figure;
subplot(3,1,1)
hold on;
stairs(1:kk, SLv.h2)
stairs(yzad,"--")
hold off
xlabel("k")
ylabel("h_2")
legend("SL","y_z_a_d")
title("Regulaotr SL, error = " + SLv.err_sum)
ylim([25, 50])

subplot(3,1,2)
stairs(1:kk,SLv.FDc)
xlabel("k")
ylabel("FD")
title("FD")
ylim([0, 37])

subplot(3,1,3)
stairs(1:kk,SLv.F1in)
xlabel("k")
ylabel("F1in")
title("F1in")
ylim([50, 90])
exportgraphics(f,'SL_zakl_zmiana_wart.pdf')

end