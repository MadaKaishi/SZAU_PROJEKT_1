function dv = stan_ciagly(t,v,F1_ust,FD,A1,C2,ap1,ap2,tau)
     
    if t>= tau
        F1 = F1_ust;
    else
        F1 = 0;
    end
    dv = zeros(2, 1);
    dv(1)= (F1 + FD - ap1 * sqrt(v(1)/A1));
    dv(2)= ap1 * sqrt(v(1)/A1) - ap2 * sqrt (sqrt(v(2)/C2));
end