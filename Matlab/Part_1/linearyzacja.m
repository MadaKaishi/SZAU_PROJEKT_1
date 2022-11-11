function dh = linearyzacja(t,h,F1_ust,FD,A1,C2,ap1,ap2,tau,h1_0,h2_0)
    if t>= tau
        F1 = F1_ust;
    else
        F1 = 78;
    end
    dh = zeros(2, 1);
    dh(1) = (F1 + FD - ap1 * (sqrt(h1_0) + (h(1) - h1_0)/(2*sqrt(h1_0))))/A1;
    dh(2) = (ap1 *((sqrt(h1_0)) + (h(1)-h1_0)/(2*sqrt(h1_0))) - ap2*(sqrt(h2_0)+ (h(2)-h2_0)/(2*sqrt(h2_0))))/(2*C2*h(2));
end