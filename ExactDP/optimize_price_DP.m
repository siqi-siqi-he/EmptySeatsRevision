function [fp1,fp2,fp3,fp2j,fp3j,r1,r2j,r3j,ofv] = optimize_price_DP(c,y,x,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,approxxyt1,approxxy1t1,approxxjy1t1,approxxjjy2t1,noExtra)
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
e0 = exp(a0);
m = optimproblem('ObjectiveSense','max');
opt = optimoptions("fmincon","Algorithm","sqp","Display","off","StepTolerance",10^-12);
p1 = optimvar('p1','LowerBound',0,'UpperBound',p1max,'Type','continuous');
p2 = optimvar('p2','LowerBound',0,'UpperBound',p2max,'Type','continuous');
p3 = optimvar('p3','LowerBound',0,'UpperBound',p3max,'Type','continuous');
p2j = optimvar('p2j',c,'LowerBound',0,'UpperBound',1,'Type','continuous');
p3j = optimvar('p3j',c,'LowerBound',0,'UpperBound',1,'Type','continuous');
m.Constraints.nonneg = p1+p2+p3<=1-p0min;
m.Constraints.sum2j = sum(p2j)==p2;
m.Constraints.sum3j = sum(p3j)==p3;
m.Constraints.excludeoccupied2 = p2j(x==0,1)==eps;
% m.Constraints.r1pos = a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))>=0;
% m.Constraints.r2pos = a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))>=0;
% m.Constraints.r3pos = a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j.*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))>=0;
con = (1-p0min)/3;
x0.p1 = max(eps,(1-p0min)/3-0.1);
x0.p2 = con+(c-sum(x))*eps;
x0.p2j = con/sum(x)*x';
x0.p2j(x==0) = eps;
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
    % m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
    % +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'))...
    % +sum(p3j.*(a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j.*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))+approxxjjy2t1'));
    m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*((a1-a0/tau1)/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1))+approxxy1t1)...
    +sum(p2j.*((a2'-a0/tau2)/b2+1/b2*(log(1-p1-p2-p3)-log(p2j))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2))+approxxjy1t1'))...
    +sum(p3j.*((a3'-a0/tau3)/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3))+approxxjjy2t1'));
else
    x0.p3 = 0;
    x0.p3j = zeros(c,1);
    m.Constraints.no3 = p3==0;
    m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
    +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'));
    % m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*((a1-a0/tau1)/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1))+approxxy1t1)...
    % +sum(p2j.*((a2'-a0/tau2)/b2+1/b2*(log(1-p1-p2-p3)-log(p2j))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2))+approxxjy1t1'));
end

[sol,ofv,ef,out] = solve(m,x0,Options=opt);
fp1 = sol.p1;
fp2 = sol.p2;
fp3 = max(0,sol.p3);
fp2j = sol.p2j;
fp3j = max(0,sol.p3j);
% r1 = a1/b1+1/(b1*tau1)*(log(1.0-fp1-fp2-fp3)-log(fp1*e0));
% r2j = a2'/b2+1/b2*(log(1-fp1-fp2-fp3)-log(fp2j*e0))+1/b2*(1-tau2)/tau2*(log(1-fp1-fp2-fp3)-log(fp2*e0));
% r3j = a3'/b3+1/b3.*(log(1-fp1-fp2-fp3)-log(fp3j*e0))+1/b3*(1-tau3)/tau3*(log(1-fp1-fp2-fp3)-log(fp3*e0));
r1 = (a1-a0/tau1)/b1+1/(b1*tau1)*(log(1.0-fp1-fp2-fp3)-log(fp1*e0));
r2j = (a2'-a0/tau2)/b2+1/b2*(log(1-fp1-fp2-fp3)-log(fp2j*e0))+1/b2*(1-tau2)/tau2*(log(1-fp1-fp2-fp3)-log(fp2*e0));
r3j = (a3'-a0/tau3)/b3+1/b3.*(log(1-fp1-fp2-fp3)-log(fp3j*e0))+1/b3*(1-tau3)/tau3*(log(1-fp1-fp2-fp3)-log(fp3*e0));
end
% 
% function [fp1,fp2,fp3,fp2j,fp3j,r1,r2j,r3j,ofv] = optimize_price_DP(c,y,x,a0,a1,a2,a3,b1,b2,b3,tau1,tau2,tau3,approxxyt1,approxxy1t1,approxxjy1t1,approxxjjy2t1,noExtra)
% echo off
% warning('off','MATLAB:nearlySingularMatrix')
% warning('off','MATLAB:SingularMatrix')
% x = x';
% 
% %bounds when all r == 0
% p0min = exp(a0)./(exp(a0)+exp(a1*tau1)+sum(exp(a2).^tau2)+sum(exp(a3).^tau3));
% p1max = exp(a1*tau1)/(exp(a0)+exp(a1*tau1));
% p2max = sum(exp(a2).^tau2)/(exp(a0)+sum(exp(a2).^tau2));
% p3max = sum(exp(a3).^tau3)/(exp(a0)+sum(exp(a3).^tau3));
% e0=exp(a0);
% xjstar = x;
% eps = 0.00000001*any([sum(x)<c-1,sum(x)~=y])+0.00000001*all([sum(x)>=c-1,sum(x)==y]);
% xjstar(1:2:c) = x(2:2:c);
% xjstar(2:2:c) = x(1:2:c);
% adj = 1*(1-(y<=1))*any(x+xjstar==2); % set double seat availability to 1, only if y>=2 and double seats are available
% %real function code
% m = optimproblem('ObjectiveSense','max');
% opt = optimoptions("fmincon","Algorithm","sqp","Display","none",'MaxFunctionEvaluations',10000,"EnableFeasibilityMode",true);
% % optimoptions("surrogateopt","Display","none");
% p1 = optimvar('p1','LowerBound',0,'UpperBound',p1max,'Type','continuous');
% p2 = optimvar('p2','LowerBound',0,'UpperBound',p2max,'Type','continuous');
% p3 = optimvar('p3','LowerBound',0,'UpperBound',p3max,'Type','continuous');
% p2j = optimvar('p2j',c,'LowerBound',0,'UpperBound',p2max*(x+eps/c),'Type','continuous');
% p3j = optimvar('p3j',c,'LowerBound',0,'UpperBound',p3max*(x+eps/c),'Type','continuous');
% if adj == 1
%     m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
%         +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'))...
%         +sum(p3j.*(a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))+approxxjjy2t1'));
%     m.Constraints.excludeoccupied3a = p3j(x==0,1)==eps;
%     m.Constraints.excludeoccupied3b = p3j(xjstar==0,1)==eps;
%     m.Constraints.posprice3 = a3'/b3+1/b3.*(log(1-p1-p2-p3)-log(p3j*e0))+1/b3*(1-tau3)/tau3*(log(1-p1-p2-p3)-log(p3*e0))>=0;
% else
%     m.Objective = (1.0-p1-p2-p3)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))+approxxy1t1)...
%         +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))+approxxjy1t1'));
%     m.Constraints.no3 = p3==0;
% end
% m.Constraints.nonneg = p1+p2+p3<=1-p0min;
% m.Constraints.sum2j = sum(p2j)==p2;
% m.Constraints.sum3j = sum(p3j)==p3;
% m.Constraints.excludeoccupied2 = p2j(x==0,1)==eps;
% m.Constraints.r1pos = a1/b1+1/(b1*tau1)*(log(1-p1-p2-p3)-log(p1*e0))>=0;
% m.Constraints.r2pos = [a2'/b2+1/b2*(log(1-p1-p2-p3)-log(p2j*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2-p3)-log(p2*e0))]'>=0;
% if noExtra == 1
%     x0.p3 = 0;
%     x0.p3j = zeros(c,1);
%     m.Constraints.no3 = p3==0;
%     r1=30;
%     r2=30;
%     x0.p1 = exp((a1-b1.*r1).*tau1)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2));
%     x02 = exp((a2-b2.*r2).*tau2)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2));
%     con = (1-p0min)/3-eps;
%     x0.p1 = max(eps,con-0.2);
%     x0.p2 = con+(c-sum(x))*eps;
%     x0.p2j = (con)/sum(x)*x';
%     x0.p2j(x==0) = eps;
%     x0.p2j(x==1) = x02(x==1);
%     x0.p2 = sum(x0.p2j);
%     m.Objective = (1.0-p1-p2)*approxxyt1+p1*(a1/b1+1/(b1*tau1)*(log(1-p1-p2)-log(p1*e0))+approxxy1t1)...
%     +sum(p2j.*(a2'/b2+1/b2*(log(1-p1-p2)-log(p2j.*e0))+1/b2*(1-tau2)/tau2*(log(1-p1-p2)-log(p2*e0))+approxxjy1t1'));
% elseif adj == 1
%     con = (1-p0min)/3-eps;
%     r1=47.5;
%     r2=47.5;
%     r3=47.5;
%     x0.p3 = con+eps*(c-sum((x+xjstar==2)));
%     x0.p3j = (con)/sum((x+xjstar==2))*(x+xjstar==2)';
%     x0.p3j(x==0) = eps;
%     x0.p3j(xjstar==0) = eps;
%     x0.p3j(x==0) = eps;
%     x0.p3j(xjstar==0) = eps;
%     x03 = exp((a3-b3.*r3).*tau3)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%     x0.p3j(x+xjstar==2) = x03(x+xjstar==2);
%     x0.p3 = sum(x0.p3j);
%     x0.p1 = exp((a1-b1.*r1)).^tau1./(exp(a0)+exp((a1-b1.*r1)).^tau1+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%     x02 = exp((a2-b2.*r2).*tau2)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%     con = (1-p0min)/3-eps;
%     x0.p1 = max(eps,con-0.2);
%     x0.p2 = con+(c-sum(x))*eps;
%     x0.p2j = (con)/sum(x)*x';
%     x0.p2j(x==0) = eps;
%     x0.p2j(x==1) = x02(x==1);
%     x0.p2 = sum(x0.p2j);
% else
%     r1=47.5;
%     r2=47.5;
%     r3=47.5;
%     x0.p3 = 0;
%     x0.p3j = zeros(c,1);
%     x0.p1 = exp((a1-b1.*r1).*tau1)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2));
%     x02 = exp((a2-b2.*r2).*tau2)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2));
%     con = (1-p0min)/3-eps;
%     x0.p1 = max(eps,con-0.2);
%     x0.p2 = con+(c-sum(x))*eps;
%     x0.p2j = (con)/sum(x)*x';
%     x0.p2j(x==0) = eps;
%     x0.p2j(x==1) = x02(x==1);
%     x0.p2 = sum(x0.p2j);
% end
% [sol,ofv,ef,output] = solve(m,x0,Options=opt);
% if ef >= 2 && isempty(output.bestfeasible)==0
%     ofv = -output.bestfeasible.fval;
%     % constrviola = output.bestfeasible.constrviolation;
%     sol = output.bestfeasible.x;
% end
% fp1 = sol.p1;
% fp2 = sol.p2;
% fp3 = sol.p3;
% fp2j = sol.p2j;
% fp3j = sol.p3j;
% r1 = a1/b1+1/(b1*tau1)*(log(1.0-fp1-fp2-fp3)-log(fp1*e0));
% r2j = a2'/b2+1/b2*(log(1-fp1-fp2-fp3)-log(fp2j*e0))+1/b2*(1-tau2)/tau2*(log(1-fp1-fp2-fp3)-log(fp2*e0));
% r3j = a3'/b3+1/b3.*(log(1-fp1-fp2-fp3)-log(fp3j*e0))+1/b3*(1-tau3)/tau3*(log(1-fp1-fp2-fp3)-log(fp3*e0));
% counter = 0;
% while any([r1;r2j;r3j]<-0.05) || ef==-2 && counter<10
% %     disp([r1;r2j;r3j]');
% %     disp(ef);
%     r1=r1*0.95;
%     r2=r2.*0.95;
%     r3=r3.*0.95;
%     x0.p1 = exp((a1-b1.*r1).*tau1)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%     x02 = exp((a2-b2.*r2).*tau2)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%     x0.p2j(x==1) = x02(x==1);
%     x0.p2 = sum(x0.p2j);
%     if adj == 1
%         x0.p3j(x==0) = eps;
%         x0.p3j(xjstar==0) = eps;
%         x03 = exp((a3-b3.*r3).*tau3)./(exp(a0)+exp((a1-b1.*r1)*tau1)+sum(exp(a2-b2.*r2).*tau2)+sum(exp(a3-b3.*r3).*tau3));
%         x0.p3j(x+xjstar==2) = x03(x+xjstar==2);
%         x0.p3 = sum(x0.p3j);
%     else
%         x0.p3 = 0;
%         x0.p3j = zeros(c,1);
%     end
%     [sol,ofv,ef,output] = solve(m,x0,Options=opt);
%     if ef >= 2 && isempty(output.bestfeasible)==0
%         ofv = -output.bestfeasible.fval;
%         % constrviola = output.bestfeasible.constrviolation;
%         sol = output.bestfeasible.x;
%     end
%     fp1 = sol.p1;
%     fp2 = sol.p2;
%     fp3 = sol.p3;
%     fp2j = sol.p2j;
%     fp3j = sol.p3j;
%     r1 = a1/b1+1/(b1*tau1)*(log(1.0-fp1-fp2-fp3)-log(fp1*e0));
%     r2j = a2'/b2+1/b2*(log(1-fp1-fp2-fp3)-log(fp2j*e0))+1/b2*(1-tau2)/tau2*(log(1-fp1-fp2-fp3)-log(fp2*e0));
%     r3j = a3'/b3+1/b3.*(log(1-fp1-fp2-fp3)-log(fp3j*e0))+1/b3*(1-tau3)/tau3*(log(1-fp1-fp2-fp3)-log(fp3*e0));
%     counter = counter+1;
% end
% if any([r1;r2j;r3j]<-0.05) 
%     disp([r1;r2j;r3j]');
% end
% end
