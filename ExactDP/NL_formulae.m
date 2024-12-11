Y = sum(c)+2;
meanprices = zeros([2^c+1,1+2*c]);
for D = 4:(2^c+1)
    x = DtoB0(D,c);
    xstar = x;
    for j = 1:2:c
        xstar([j j+1]) = x([j+1 j]);
    end
    x2 = [1*(sum(x)>0),x,x.*xstar];
    YD = sum(x)+2;
    if YD == 3
        val = (squeeze(policy2(D,3:YD,:,2:17)).*x2')';
    else
        val = squeeze(policy2(D,3:YD,:,2:17)).*x2;
    end 
    mask = (val<inf) & (val>0);
    meanprices(D,:) = median(val.*mask,[1,3],'omitnan');
end
meanval = zeros(2*c+1,1);
for seat = 1:(2*c+1) %seat does not mean seat, but seat product combination
    meanpricesseat = meanprices(:,seat);
    meanval(seat) = median(meanpricesseat(meanpricesseat~=0),'omitnan');
end

% %when noExtra1 loaded
% file_name = sprintf('meanprices_a0_%d_noExtra1.mat',a0);
% save(file_name,"meanval");
% file_name = sprintf('meanprices_a0_%d_noExtra1.txt',a0);
% writematrix(meanval, file_name);

%when noExtra0 loaded
file_name = sprintf('meanprices_a0_%d_median.mat',a0);
save(file_name,"meanval");
file_name = sprintf('meanprices_a0_%d_median.txt',a0);
writematrix(meanval, file_name);
