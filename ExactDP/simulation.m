function [rev_nt,policy,state,tickets_n,seats2_n,seats3_n] = simulation(nsim,dec,T,c,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra)
rev_nt = zeros(nsim,T);
policy = zeros(nsim,T,1+2*c);
state = zeros(nsim,c+1);
tickets_n = zeros(nsim,4);
seats2_n = zeros(nsim,c);
seats3_n = zeros(nsim,c);
parfor n = 1:nsim 
    x = ones(c,1);
    y = sum(x);
    cTickets = zeros(4,1);
    cSeats2 = zeros(c,1);
    cSeats3 = zeros(c,1);
    for t = 1:T
        t2 = T+2-t;
        if sum(x)==0 || y==0
            % break
        else
            [fp1,~,~,fp2j,fp3j,r1,r2j,r3j,~] = opt_DP_returnProbs(c,x,y,t2,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,V,noExtra);
            policy(n,t,:) = [r1;r2j;r3j];
            [x,y,rev,cTickets,cSeats2,cSeats3] = calc_decision(fp1,fp2j,fp3j,r1,r2j,r3j,dec,x,y,t,n,cTickets,cSeats2,cSeats3);
            rev_nt(n,t) = max(0,rev);
        end
    end
    tickets_n(n,:) = cTickets;
    seats2_n(n,:) = cSeats2;
    seats3_n(n,:) = cSeats3;
end
end