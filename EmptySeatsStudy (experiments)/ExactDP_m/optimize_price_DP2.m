function [fp1,fp2,fp3,fp2j,fp3j,r1,r2j,r3j,ofv] = optimize_price_DP2(c,y,x,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,approxxyt1,approxxy1t1,approxxjy1t1,approxxjjy2t1,noExtra)
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:SingularMatrix')
%bounds when all r == 0
p0min = exp(a0)./(exp(a0)+exp(a1*tau1)+sum(exp(a2).^tau2)+sum(exp(a3).^tau3));
p1max = exp(a1*tau1)/(exp(a0)+exp(a1*tau1));
p2max = sum(exp(a2).^tau2)/(exp(a0)+sum(exp(a2).^tau2));
p3max = sum(exp(a3).^tau3)/(exp(a0)+sum(exp(a3).^tau3));
%data preparation
xjstar = x;
eps = 0.00000001;
for j = 1:2:c
    xjstar([j j+1]) = x([j+1 j]);
end
%adj only available when adjacent seats unoccupied and y>1
adj = any(x+xjstar==2)*(1-(y<=1)); 
%model
m = optimproblem('ObjectiveSense','max');
opt = optimoptions("fmincon","Algorithm","interior-point","Display","off","StepTolerance",10^-30,"MaxFunctionEvaluations",9*10^5,"MaxIterations",9*10^4);
p1 = optimvar('p1','LowerBound',0,'UpperBound',p1max,'Type','continuous');
p2 = optimvar('p2','LowerBound',0,'UpperBound',p2max,'Type','continuous');
p3 = optimvar('p3','LowerBound',0,'UpperBound',p3max,'Type','continuous');
p2j = optimvar('p2j',c,'LowerBound',0,'UpperBound',1,'Type','continuous');
p3j = optimvar('p3j',c,'LowerBound',0,'UpperBound',1,'Type','continuous');
r2 = optimvar('r2',1,'LowerBound',0,'Type','continuous');
r3 = optimvar('r3',1,'LowerBound',0,'Type','continuous');
m.Constraints.nonneg = p1+p2+p3<=1-p0min;
m.Constraints.sum2j = sum(p2j)==p2;
m.Constraints.sum3j = sum(p3j)==p3;
m.Constraints.excludeoccupied2 = p2j(x==0,1)==eps;
con = (1-p0min)/3;
x0.p1 = max(eps,(1-p0min)/3-0.1);
x0.p2 = con+(c-sum(x))*eps;
x0.p2j = con/sum(x)*x';
x0.p2j(x==0) = eps;
e0 = exp(a0);
if noExtra == 1
    x0.p3 = 0;
    x0.p3j = zeros(c,1);
    m.Constraints.no3 = p3==0;
    m.Objective = (1.0-p1-p2)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2)-log(p1*e0))+approxxy1t1)...
    +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2)-log(p2*e0))+approxxjy1t1'));
elseif adj == 1
    x0.p3 = con+eps*(c-sum((x+xjstar==2)));
    x0.p3j = con/sum((x+xjstar==2))*(x+xjstar==2)';
    x0.p3j(x==0) = eps;
    x0.p3j(xjstar==0) = eps;
    m.Constraints.excludeoccupied3a = p3j(x==0,1)==eps;
    m.Constraints.excludeoccupied3b = p3j(xjstar==0,1)==eps;
    m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
    +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'))...
    +sum(p3j.*(a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j.*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))+approxxjjy2t1'));
else
    x0.p3 = 0;
    x0.p3j = zeros(c,1);
    m.Constraints.no3 = p3==0;
    m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
    +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'));
end

%set one price per product
m.Constraints.prices2 = (a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))).*(x)'==r2.*(x)';
m.Constraints.prices3 =  (a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j.*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))).*(xjstar+x==2)'==r3.*(xjstar+x==2)';
x0.r2 = 20;
x0.r3 = 20;
x0.p2j = exp((a2-b2*x0.r2).*tau2)./(exp(a0)+exp((a1-b1*10)*tau1)+sum(exp(a2-b2*x0.r2).*tau2)+sum(exp(a3-b3*x0.r3).*tau3));
x0.p3j = exp((a3-b3*x0.r3).*tau3)./(exp(a0)+exp((a1-b1*10)*tau1)+sum(exp(a2-b2*x0.r2).*tau2)+sum(exp(a3-b3*x0.r3).*tau3));
x0.p2 = sum(x0.p2j);
x0.p3 = sum(x0.p3j);

[sol,ofv,ef,out] = solve(m,x0,Options=opt);
fp1 = sol.p1;
fp2 = sol.p2;
fp3 = max(0,sol.p3);
fp2j = sol.p2j;
fp3j = max(0,sol.p3j);
r1 = a1/b1+1/(b1*tau1)*(log(1.0-fp1-fp2-fp3)-log(fp1*e0));
r2j = a2'/b2+1/b2*(log(1-fp1-fp2-fp3)-log(fp2j*e0))+1/b2*(1-tau2)/tau2*(log(1-fp1-fp2-fp3)-log(fp2*e0));
r3j = a3'/b3+1/b3.*(log(1-fp1-fp2-fp3)-log(fp3j*e0))+1/b3*(1-tau3)/tau3*(log(1-fp1-fp2-fp3)-log(fp3*e0));
end