%% Plot
genes_ranking = [Breast7390(1:200),Colon17536(1:200),Liver10141(1:200),Lung11969(1:200),Osteosarcoma21257(1:200),Ovarian26712_netcox(1:200),Renal29609(1:200)];
TopRankedGenes = TopRanked(genes_ranking);
N = size(TopRankedGenes,2);

figure
hold all;
ColorSet = varycolor(N); 
set(gca, 'FontSize', 12)
set(gca, 'ColorOrder', ColorSet); 
box;
for i =1:N
    plot(TopRankedGenes(:,i),'LineWidth', 2);
end
xlabel('Number of Top-ranked Genes')  
ylabel('Percentage of Overlapped Genes')
legend('Breast','Colon', 'Liver', 'Lung', 'Osteosarcoma', 'Ovarian', 'Renal')