c = 8;
x = ones(c);
pr0 = 1;
pr1 = 1;
pr2j = 1;
pr3j = 1;
for t = T+1:-1:1
    pol = policy2();
    r1 = pol(1);
    r2 = pol(2:c+1);
    r3 = pol(c+2:end);
    p1 = exp((a1-b1.*r1).*tau1)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
    p2j = exp((a2-b2.*r2).*tau2)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
    p3j = exp((a3-b3.*r3).*tau3)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
    pr1 = pr1.*p1;
    pr2j = pr2j.*p2j;
    pr3j = pr3j.*p3j;

end
for y = 3:10
    polW3 = policyWITH3(1:257,y,2,2:end);
    polWO3 = policyWITHOUT3(1:257,y,2,2:end);
    pW3 = polW3(polW3~=0);
    pWO3 = polWO3(polWO3~=0);
    mean((pW3-pWO3)./pWO3,'all')
end
for y = 3:10
    polW3 = policyWITH3(1:257,y,1,2:end);
    polWO3 = policyWITHOUT3(1:257,y,1,2:end);
    pW3 = polW3(polW3~=0);
    pWO3 = polWO3(polWO3~=0);
    mean((pW3-pWO3)./pWO3,'all')
end 