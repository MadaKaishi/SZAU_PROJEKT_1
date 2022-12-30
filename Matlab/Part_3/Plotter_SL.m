clear all

porown = false;

addpath("../Part_3")

generate = false;


if generate
clear all


DMC_new_traj
save("Do_porownania\DMC.mat", "h2","F1in","error");
clear all

FDMC_new_traj
save("Do_porownania\FDMC.mat", "h2","F1in","err_sum");
clear all

end

ks = 1500;
kp = 1402;
kk = 5000;

yzad(1:ks)=38.44; 
yzad(ks:5000)=80; 

draws = true;
if porown
DMCv = load("Do_porownania/DMC.mat");
FDMCv = load("Do_porownania/FDMC.mat");
SLv = load("Pliki_mat_SL_norm/SL_norm.mat");

close all

if draws
%Plot wyjście
figure;
stairs(DMCv.h2)
hold on;
stairs(FDMCv.h2)
stairs(SLv.h2)
stairs(yzad,"--");
hold off;
xlabel('k'); ylabel("h");
legend("DMC","FDMC","SL","h_2_z_a_d")
title("Regulaotr DMC, error = " + DMCv.error + newline +"Regulator FDMC, error = " + FDMCv.err_sum + newline +"Regulator SL, error = " + SLv.err_sum)
exportgraphics(gca,'SL_porown.pdf')

%Plot sterowanie
figure;
stairs(DMCv.F1in)
hold on
stairs(FDMCv.F1in)
stairs(SLv.F1in)
legend("DMC","FDMC","SL")
xlabel('k'); ylabel("F_1_i_n");
title("Sterowanie regulatorów")
exportgraphics(gca,'SL_porown_ster.pdf')

end


else

SLv = load("Pliki_mat_SL_norm/SL_norm.mat");
%Plot wyjście
figure;
stairs(SLv.h2)
hold on
stairs(yzad,"--");
hold off;
xlabel('k'); ylabel("h");
legend("SL","h_2_z_a_d")
title("Regulator SL, error = " + SLv.err_sum)
exportgraphics(gca,'SL_zmana_wart.pdf')

%Plot sterowanie
figure;
stairs(SLv.F1in)
legend("SL")
xlabel('k'); ylabel("F_1_i_n");
title("Sterowanie regultorów")
exportgraphics(gca,'SL_zmana_wart_ster.pdf')

end