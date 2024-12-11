function [a0,a1,a2,a3,b1,b2,b3,tau0,tau1,tau2,tau3] = pricing_homo_seats(c)
    a0 = 0;
    a1 = 0.2;
    a2 = 0.4 * ones(1,c);
    a3 = 0.6 * ones(1,c);
    b1 = 0.2;
    b2 = 0.4;
    b3 = 0.6;
    tau0 = 1;
    tau1 = 1;
    tau2 = 0.4;
    tau3 = 0.6;
end