function [fp1,fp2,fp3,fp2j,fp3j,r1,r2j,r3j,ofv] = opt_DP_returnProbs(c,x,y,t,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra)
approxxjy1t1 = zeros(1,c);
approxxjjy2t1 = zeros(1,c);
x2j = x'-eye(c);
x3j = x'-eye(c);   
y2 = y+2; %Index does not equal the number of possible sales
approxxyt1 = V(BtoD0(x),y2,t-1);
approxxy1t1 = -9999*(y-1<0)+V(BtoD0(x),y2-1,t-1)*(y-1>=0);
for j = 1:c
    x3j(j,j-(-1)^(j)) = x3j(j,j-(-1)^(j))-1;
    approxxjy1t1(j) = -9999*(y-1<0)+V(BtoD0(x2j(j,:)),y2-1,t-1)*(y-1>=0);
    approxxjjy2t1(j) = -9999*(y-2<0)+V(BtoD0(x3j(j,:)),y2-2,t-1)*(y-2>=0);
end
% Solve in MATLAB
[fp1,fp2,fp3,fp2j,fp3j,r1,r2j,r3j,ofv] = optimize_price_DP(c,y,x,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,approxxyt1,approxxy1t1,approxxjy1t1,approxxjjy2t1,noExtra);
end