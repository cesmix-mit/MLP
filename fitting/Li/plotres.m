filename = "numdesc.txt";
fileID = fopen(filename,'r');
num = fscanf(fileID, '%f');
fclose(fileID);

num = reshape(num,[3, 6])';

filename = "testerror.txt";
fileID = fopen(filename,'r');
err = fscanf(fileID, '%f');
fclose(fileID);
err = reshape(err,[6, 6])';

err(:,1:6) = 1e3*err(:,1:6);
figure(1); clf; 
semilogy(num(:,1), err(:,1),'-.o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on
semilogy(num(:,2), err(:,2), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), err(:,3), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca, 'YTick', [0.5 1 2 4 8 16])
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 18, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 145 0.5 16]);
text(70,8.0,"Li","FontSize",28,"FontWeight","bold");
print -dpng Li_test_energy_error.png

figure(2); clf; 
semilogy(num(:,1), err(:,4), '-.o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on;
semilogy(num(:,2), err(:,5), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), err(:,6), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca, 'YTick', [10 20 40 80])
% xlabel("Number of descriptors", 'FontSize', 18);
% ylabel("MAE force error (meV/Å)", 'FontSize', 18);
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 18, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 145 10 100]);
text(70,0.065*1e3,"Li","FontSize",28,"FontWeight","bold");
print -dpng Li_test_force_error.png



