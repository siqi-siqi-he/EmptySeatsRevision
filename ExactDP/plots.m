% Data from the table
a0 = [-1, -0.5, 0, 0.5, 1];
Relative_Exp = [100.61, 101.52, 103.45, 106.88, 112.03];
SLF_with_P3 = [0.85349179,0.79692923,0.72400082,0.63967122,0.54995664]*100;  
SLF_without_P3 = [0.85953338,0.80497636,0.72662882,0.62247146,0.50103395]*100;  
 

% Plot lines with dots
figure1 = figure('Color',[1 1 1]);
plot(a0, Relative_Exp, '-o', a0, SLF_with_P3, '-^', a0, SLF_without_P3, '-s');

% Add labels and title
xlabel('Demand');
ylabel('Percentage');
% title('Increase of Relative Expected Revenue by selling Extra Seats over Demand Ratio (HOMOG)');
xticks(a0);
xticklabels({'Strong, a0=-1','','','','Weak, a0=1'})

% Add legend
legend('Relative Expected Revenue', 'SLF with Extra Seats', 'SLF without Extra Seats', 'Location', 'Best');

% Adjust axes limits if necessary
xlim([min(a0)-0.3, max(a0)+0.3]);
ylim([min([Relative_Exp, SLF_with_P3, SLF_without_P3])-1, max([Relative_Exp, SLF_with_P3, SLF_without_P3])+1]);
f = gcf;
saveas(f,'Fig5','epsc');





% Generate sample data for histograms
D=3;
DtoB0(D,8)
data1 = policy2(D,:,1,:);  % Sample data for histogram 1
data2 = policy2(D,:,2,:);  % Sample data for histogram 2
data3 = policy2(D,:,10,:);  % Sample data for histogram 3
data1 = data1(data1~=inf);
data2 = data2(data2~=inf);
data3 = data3(data3~=inf);

% Create a new figure
figure;

% Plot histogram 1
subplot(1, 3, 1);  % Create subplot in a 1x3 grid, this is the first subplot
histogram(data1);
title('Product 1');
xlabel('Prices');
ylabel('Frequency');

% Plot histogram 2
subplot(1, 3, 2);  % Create subplot in a 1x3 grid, this is the second subplot
histogram(data2);
title('Product 2');
xlabel('Prices');
ylabel('Frequency');

% Plot histogram 3
subplot(1, 3, 3);  % Create subplot in a 1x3 grid, this is the third subplot
histogram(data3);
title('Product 3');
xlabel('Prices');
ylabel('Frequency');