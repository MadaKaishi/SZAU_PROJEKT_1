function h2 = fun_ob(var)

%Constants
A1 = 505;
C2 = 0.65;
ap1 = 23; %alfa_1
ap2 = 15; %alfa_2

%Punkt pracy
tau = 120;
h2_0 = var(1);
h1_0 = var(2);
F1 = var(3);
FD = 15;

t_sym = 5000; %czas symulacji
T = 1; %krok

%warunki_początkowe
kp = tau/T + 2;
kk = t_sym/T;
h1(1:kp) = h1_0;
h2(1:kp) = h2_0;
F1in(1:kk) = F1;
FDc(1:kk) = FD;


for k = kp:kk
      h1(k) = h1(k-1) + T*(F1in(k-1-(tau/T)) + FDc(k-1) - ap1 * (sqrt(h1_0) + (h1(k-1) - h1_0)/(2*sqrt(h1_0))))/A1;
      h2(k) = h2(k-1) + T*((ap1 * sqrt(h1_0) - ap2 * sqrt(h2_0))/(2*C2*h2_0)   +    ap1*(h1(k-1) - h1_0)/(4*C2*(sqrt(h1_0)))    -   ap2*(h2(k-1) - h2_0)/(4*C2*(nthroot((h2_0 * h2_0),3))));
end

