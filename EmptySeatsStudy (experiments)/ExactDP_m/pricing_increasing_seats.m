function [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_increasing_seats(c)
    add = linspace(0,0.01*c-0.01,c);addhelp = (add(2:2:c)+add(1:2:c))./2;
    add2 = repmat(addhelp',1,2)';
    add2 = add2(:)';
    a0 = 0;
    a1 = 0.2;
    a2 = 0.4 * ones(1,c) + add;
    a3 = 0.6 * ones(1,c) + add2;
    b1 = 0.2;
    b2 = 0.4;
    b3 = 0.6;
    tau0 = 1;
    tau1 = 1;
    tau2 = 0.4;
    tau3 = 0.6;
end