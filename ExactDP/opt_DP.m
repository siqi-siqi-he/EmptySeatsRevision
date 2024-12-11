function [ofv,r1,r2j,r3j,SLFx] = opt_DP(V,c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,noExtra,SLF)
    approxxyt1 = V(BtoD0(x),y,t-1);
    approxxy1t1 = V(BtoD0(x),y-1,t-1);
    approxxjy1t1 = zeros(1,c);
    approxxjjy2t1 = zeros(1,c);
    x2j = x-eye(c);
    x3j = x-eye(c);
    for j = 1:c
        x3j(j,j-(-1)^(j)) = x3j(j,j-(-1)^(j))-1; 
        approxxjy1t1(j) = V(BtoD0(x2j(j,:)),y-1,t-1);
        approxxjjy2t1(j) = V(BtoD0(x3j(j,:)),y-2,t-1);
    end
    y2 = y - 2; %Index does not equal the number of possible sales
    [p1,p2,p3,p2j,p3j,r1,r2j,r3j,ofv] = optimize_price_DP(c,y2,x,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,approxxyt1,approxxy1t1,approxxjy1t1,approxxjjy2t1,noExtra);
    SLFx = (1-p1-p2-p3)*SLF(BtoD0(x),y,t-1)+p1+p2+p3+p1*SLF(BtoD0(x),y-1,t-1);
    for j = 1:c
        SLFx = SLFx+p2j(j).*SLF(BtoD0(x2j(j,:)),y-1,t-1)+p3j(j).*SLF(BtoD0(x3j(j,:)),y-2,t-1);
    end
end