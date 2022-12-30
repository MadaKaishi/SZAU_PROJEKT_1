clear all



addpath("../Part_3")

generate = true;

if generate
    clear all
    
    DMC_new_traj
    save("Do_porownania\DMC.mat", "h2","F1in","error");
    clear all
    
    FDMC_new_traj
    save("Do_porownania\FDMC.mat", "h2","F1in","err_sum");
    clear all
    close all

end

ks = 1500;
kp = 1402;
kk = 5000;

yzad(1:ks)=38.44; 
yzad(ks:5000)=80; 

porown = true;

if porown
    DMCv = load("Do_porownania/DMC.mat");
    FDMCv = load("Do_porownania/FDMC.mat");
    SLv = load("../Part_3/Pliki_mat_SL_norm/SL_norm.mat");
    NPLv = load("NPL_norm.mat");

    close all
    
    %Plot wyjście
    figure;
    stairs(DMCv.h2)
    hold on;
    stairs(FDMCv.h2)
    stairs(SLv.h2)
    stairs(NPLv.h2)
    stairs(yzad,"--");
    hold off;
    xlabel('k'); ylabel("h");
    legend("DMC","FDMC","SL","h_2_z_a_d")
    title("Regulaotr DMC, error = " + DMCv.error + newline +"Regulator FDMC, error = " + FDMCv.err_sum + newline +"Regulator SL, error = " + SLv.err_sum + newline +"Regulator NPL, error = " + NPLv.err_sum)
    exportgraphics(gca,'NPL_porown.pdf')
    
    %Plot sterowanie
    figure;
    stairs(DMCv.F1in)
    hold on
    stairs(FDMCv.F1in)
    stairs(SLv.F1in)
    stairs(NPLv.F1in)
    legend("DMC","FDMC","SL","NPL")
    xlabel('k'); ylabel("F_1_i_n");
    title("Sterowanie regulatorów")
    exportgraphics(gca,'NPL_porown_ster.pdf')

else

    NPLv = load("NPL_norm.mat");
    %Plot wyjście
    figure;
    stairs(NPLv.h2)
    hold on
    stairs(yzad,"--");
    hold off;
    xlabel('k'); ylabel("h");
    legend("NPL","h_2_z_a_d")
    title("Regulator NPL, error = " + NPLv.err_sum)
    exportgraphics(gca,'NPL_zmana_wart.pdf')
    
    %Plot sterowanie
    figure;
    stairs(NPLv.F1in)
    legend("NPL")
    xlabel('k'); ylabel("F_1_i_n");
    title("Sterowanie")
    exportgraphics(gca,'NPL_zmana_wart_ster.pdf')

end