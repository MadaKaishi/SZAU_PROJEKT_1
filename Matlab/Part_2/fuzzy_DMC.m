
%% Nastawy regulatora
D = 1500;
DZ = 200;
N = 125;
Nu = 2;
lambda = 15;
a = 10;
c = [8.5 17.5 24.5 33.5];
draw = true;

    % zmienne i macierze regulatora
    D = min(D, 400);
    DZ = min(DZ, 400);
    N = min(min(N,D),DZ);
    Nu = min(Nu, N);

    il = 5;
    lambda = lambda*ones(1,il);
    ku = zeros(il,D-1);
    kz = zeros(il,DZ);
    ke = zeros(1,il);
    ymin = 4.5;
    ymax = 34.5;
%     if length(a) == 0
%         [a, c, hr0] = z2_modelroz(il, false);
%     else
        hr0(1) = (c(1)+ymin)/2-1;
        hr0(il) = min((ymax+c(il-1))/2+1, ymax);
        if il > 2
            hr0(2:il-1) = (c(2:il-1)+c(1:il-2))./2;
        end
    end
%     a = 10;
%     c = [8.5000   17.5000   24.5000   33.5000];
    
    if draw
        figure
        hold on
        for r = 1:il
            if r == 1
                plot(ymin:0.1:ymax,1-1./(1+exp(-a*([ymin:0.1:ymax]-c(1)))))
            elseif r == il
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(il-1)))))
            else
                plot(ymin:0.1:ymax,1./(1+exp(-a*([ymin:0.1:ymax]-c(r-1))))-1./(1+exp(-a*([ymin:0.1:ymax]-c(r)))))
            end
        end
        plot(hr0, ones(1,il), 'ko')
        xlabel('h2')
        ylabel('przynale¿noœæ')
        title('Funkcje przynale¿noœci regulatorów lokalnych w regulacji rozmytej wzglêdem wartoœci wyjœcia')
    end
    
    for r = 1:il
        [s, z] = z1_step(hr0(r), false);

        M=zeros(N,Nu);
        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j)=s(i-j+1);
              end
           end
        end

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

        MZP=zeros(N,DZ);
        for i=1:N
            MZP(i,1) = z(i);
           for j=2:DZ
              if i+j-1<=DZ
                 MZP(i,j)=z(i+j-1)-z(j);
              else
                 MZP(i,j)=z(DZ)-z(j);
              end      
           end
        end

        I=eye(Nu);
        K=((M'*M+lambda(r)*I)^-1)*M';
        ku(r,:)=K(1,:)*MP;
        kz(r,:)=K(1,:)*MZP;
        ke(r)=sum(K(1,:));
    end
    deltaup=zeros(1,D-1);
    deltazp=zeros(1,DZ-1);
    deltauk = zeros(il,1);
    Umax = 90;
    Umin = 20;
    dumax = 1;
    w = zeros(1,il);

    % dane
    C1 = 0.95;
    C2 = 0.95;
    a1 = 16;
    a2 = 16;
    n = 12100;
    tau = 50;
    U0 = 54;
    D0 = 10;
    Y0 = 16;
    V1 = C1*Y0^2;
    h1 = Y0;
    V2 = V1;
    start = 100;
    dY = [start 34; start+3000 4.5; start+6000 24; start+9000 14];
    dZ = [start+1500 15; start+4500 5; start+7500 12.5; start+10500 7.5];

    U = U0*ones(1,n);
    Dist = D0*ones(1,n);
    Y = Y0*ones(1,n);
    Yz = Y;
    for i = 1:length(dY)
        Yz(dY(i,1):n) = dY(i,2);
        Dist(dZ(i,1):n) = dZ(i,2);
    end
    e = zeros(1,n);

    hold on
    for k = start:n
        % symulacja
        V1 = V1 + U(k-1-tau) + Dist(k-1) - a1*h1^0.5;
        V2 = V2 + a1*h1^0.5 - a2*Y(k-1)^0.5;
        h1 = (V1/C1)^0.5;
        Y(k) = (V2/C2)^0.5;
        % uchyb
        e(k) = Yz(k) - Y(k);

        %uwzglêdnianie zak³ócenia
        for i = DZ:-1:2
           deltazp(i) = deltazp(i-1);
        end
        deltazp(1) = Dist(k) - Dist(k-1);

        % Prawo regulacji
        for i = 1:il
            deltauk(i) = ke(i)*e(k)-ku(i,:)*deltaup'-kz(i,:)*deltazp';
            if i == 1
                w(i) = 1-1/(1+exp(-a*(Y(k)-c(1))));
            elseif i == il
                w(i) = 1/(1+exp(-a*(Y(k)-c(il-1))));
            else
                w(i) = 1/(1+exp(-a*(Y(k)-c(i-1)))) - 1/(1+exp(-a*(Y(k)-c(i))));
            end
        end
        DELTAuk = w*deltauk/sum(w);
        DELTAuk = max(min(DELTAuk, dumax),-dumax);

        for i = D-1:-1:2
          deltaup(i) = deltaup(i-1);
        end
        deltaup(1) = DELTAuk;
        U(k) = U(k-1)+deltaup(1);
        U(k) = max(min(U(k),Umax),Umin);
    end

    if draw
        figure
        subplot(3,1,1)
        plot(1:n, Yz, 'r')
        xlabel('t[s]')
        ylabel('h2[cm]')
        hold on
        plot(1:n, Y, 'b')
        subplot(3,1,2)
        plot(1:n, U)
        xlabel('t[s]')
        ylabel('F1in[cm3/s]')
        subplot(3,1,3)
        plot(1:n, Dist)
        xlabel('t[s]')
        ylabel('FD[cm3/s]')
    end
    E = sum(e.^2)/n;