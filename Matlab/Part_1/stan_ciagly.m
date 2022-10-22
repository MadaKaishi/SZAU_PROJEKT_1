function dh = stan_ciagly(t,h,F1_ust,FD,A1,C2,ap1,ap2,tau)
     
    if t>= tau
        F1 = F1_ust;
    else
        F1 = 0;
    end
    dh = zeros(2, 1);
    dh(1)= (F1 + FD - ap1*(sqrt(h(1))))/A1;
    %dh(2)= (ap1*(sqrt(h(1))) - ap2*nthroot(h(2),4))/C2;
    dh(2)= (ap1*sqrt(h(1)) - ap2*sqrt(h(2)))/C2;
end