%% Plot survival cruve

load('groups.mat')

% Breast cancer
h1 = figure
logrank(HighRisk1456,LowRisk1456)

hold all
box
title('Breast cancer','FontName','Arial','FontSize',30);
ylabel('Estimated survival function','FontName','Arial','FontSize',30); 
xlabel('Month','FontName','Arial','FontSize',30); 
%set(gca,'FontSize',30);