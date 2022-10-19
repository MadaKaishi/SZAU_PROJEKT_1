function dy = odefun(t,y)
%Constants
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy
F1 = 78;
FD = 15;

dy = zeros(2,1);
dy(1) = F1 + FD - ap1 * sqrt(dy(1)/A1);
dy(2) = ap1 * sqrt(dy(1)/A1) - ap2 * sqrt(dy(2)/C2);
end