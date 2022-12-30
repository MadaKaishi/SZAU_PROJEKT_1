
generate = false;


if generate
clear all
addpath("../Part_1")
addpath("../Part_1/Odp_skok")
addpath("../Part_2")
addpath("../Part_2/Do_porownania_zakl")

DMC_zmiana_zakl
save("Do_porownania_zakl/DMC.mat", "h2","F1in","error");
clear all

FDMC_v2_zakl
save("Do_porownania_zakl/FDMC.mat", "h2","F1in","err_sum");
clear all

end

ks = 1500;
kp = 1402;
kk = 15000;

%Skok wartosci zadanej:
yzad(1:5000)=38.44; 
FDc(1:ks) = FD;
FDc(ks:3250) = FD+7.5;
FDc(3250:5000) = FD-7.5;


draws = true;

DMCv = load("Do_porownania_zakl/DMC.mat");
FDMCv = load("Do_porownania_zakl/FDMC.mat");

if draws
%Plot wyjście
f = figure;
subplot(3,1,1)
stairs(1:kk, DMCv.h2)
hold on;
stairs(1:kk, FDMCv.h2)
hold off
xlabel("k")
ylabel("h_2")
legend("DMC","FDMC",Location="north",Orientation="horizontal")
title("Regulaotr DMC, error = " + DMCv.error + newline +"Regulator FDMC, error = " + FDMCv.err_sum )
ylim([25, 50])

subplot(3,1,2)
stairs(1:kk,FDc)
xlabel("k")
ylabel("FD")
title("FD")
ylim([0, 37])

subplot(3,1,3)
stairs(1:kk,DMCv.F1in)
hold on
stairs(1:kk,FDMCv.F1in)
hold off
xlabel("k")
ylabel("F1in")
title("F1in")
legend("DMC","FDMC",Location="south",Orientation="horizontal")
ylim([50, 90])
exportgraphics(f,'DMC_zakl_por.pdf')
end