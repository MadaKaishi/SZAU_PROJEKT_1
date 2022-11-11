clear all
clc

%Parametry programu
il_fun = 3;
draw = true;
sa = false;

%Zmienne zadaniowe
h_min = 0;
h_max = 90;
h = (h_min:1:h_max)';

nach = 3; %nachylenie funkcji 

d = (h_max-h_min)/il_fun; %szerokości funkcji przynależnośći
c = h_min+d:d:h_max-d; %punkty przegięcia

%Wybranie punktu pracy
hr0 = ones(1,il_fun);
hr0(1) = d/2;
hr0(il_fun) = min((h_max+c(il_fun-1))/2+1, h_max);
    if il_fun > 2
        hr0(2:il_fun-1) = (c(2:il_fun-1)+c(1:il_fun-2))./2;
    end


if draw
    figure
    hold on
    %Plotter funkcji przynaleznosci
    for i = 1:il_fun
        if i == 1
            plot(trapmf(h,[0 0 c(1)-nach/2 c(1)+ nach/2]));
        elseif i == il_fun
            plot(trapmf(h,[c(il_fun-1)-nach/2 c(il_fun-1)+nach/2 h_max h_max]));
        else
            plot(trapmf(h,[c(i-1)-nach/2 c(i-1)+ nach/2 c(i)-nach/2 c(i)+ nach/2]));
        end
    end
    xlim([0 90])
    xlabel("h_2"); ylabel("Funkcja przynależności");
    title(sprintf("Funkcja przynaleznosci dla %i zbiorów rozmytych",il_fun))
    if sa
        print(sprintf('funkcja_przynelznosci_%i.png',il_fun),'-dpng','-r400')
    end
end

zmienna = trapmf(h,[0 0 c(1)-nach/2 c(1)+ nach/2])

