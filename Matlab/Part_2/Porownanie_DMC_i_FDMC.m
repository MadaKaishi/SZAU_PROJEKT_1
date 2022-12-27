
generate = false;


if generate
clear all
addpath("../Part_1")
addpath("../Part_1/Odp_skok")
addpath("../Part_2")
addpath("../Part_2/Do_porownania")

DMC
save("Do_porownania\DMC.mat", "h2","F1in","error");
clear all

FDMC_v2
save("Do_porownania\FDMC.mat", "h2","F1in","err_sum");
clear all

end

ks = 1500;
kp = 1402;
kk = 20000;

yzad(1:ks)=38.44; 
yzad(ks:5000)=30;
yzad(5000:10000)=80;
yzad(10000:15000)=20;
yzad(15000:20000)=40;

draws = true;

DMCv = load("Do_porownania/DMC.mat");
FDMCv = load("Do_porownania/FDMC.mat");

if draws
iteracja = 0:1:kk-1;
%Plot wyjście
figure;
stairs(iteracja, DMCv.h2)
hold on;
stairs(iteracja, FDMCv.h2)
stairs(iteracja, yzad,"--");
hold off;
xlabel('k'); ylabel("h");
legend("DMC","FDMC","h_2_z_a_d")
title("Regulaotr DMC, error = " + DMCv.error + newline +"Regulator FDMC, error = " + FDMCv.err_sum )
% exportgraphics(gca,'DMC_rozm_zmiana_wart.pdf')

%Plot sterowanie
figure;
stairs(iteracja, DMCv.F1in)
hold on
stairs(iteracja, FDMCv.F1in)
legend("DMC","FDMC")
xlabel('k'); ylabel("F_1_i_n");
title("Sterowanie regulaotrów")
% exportgraphics(gca,'DMC_rozm_zmiana_ster.pdf')
end