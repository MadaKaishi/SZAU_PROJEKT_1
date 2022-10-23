function dv = linearyzacja(t,v,F1_ust,FD,A1,C2,ap1,ap2,tau,V1_0,V2_0)
    if t>= tau
        F1 = F1_ust;
    else
        F1 = 78;
    end
    dv = zeros(2, 1);
    dv(1) = (F1 + FD - ap1 * (sqrt(V1_0/A1) + 1/(2*A1*sqrt(V1_0/A1)) * (v(1) - V1_0)));
    dv(2) =  ap1 * (sqrt(V1_0/A1) + 1/(2*A1*sqrt(V1_0/A1)) * (v(1) - V1_0)) - ap2 * (nthroot(V2_0/C2,4) + 1/(4*C2*((V2_0/C2)^(3/4))) * (v(2) - V2_0));
end