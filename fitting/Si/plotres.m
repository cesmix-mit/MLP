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
set(gca, 'YTick', [4 8 16 32 64])
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 18, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 145 4 64]);
text(70,37,"Si","FontSize",28,"FontWeight","bold");
print -dpng Si_test_energy_error.png

figure(2); clf; 
semilogy(num(:,1), err(:,4), '-.o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on;
semilogy(num(:,2), err(:,5), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), err(:,6), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca, 'YTick', [90 130 180 250 360])
% xlabel("Number of descriptors", 'FontSize', 18);
% ylabel("MAE force error (meV/Å)", 'FontSize', 18);
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 18, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 145 90 360]);
text(70,0.275*1e3,"Si","FontSize",28,"FontWeight","bold");
print -dpng Si_test_force_error.png


filename = "testerror.txt";
fileID = fopen(filename,'r');
err = fscanf(fileID, '%f');
fclose(fileID);
err1 = reshape(err,[6, 6])';

filename = "trainerror.txt";
fileID = fopen(filename,'r');
err = fscanf(fileID, '%f');
fclose(fileID);
err2 = reshape(err,[6, 6])';

e1 = 1e3*err1(:,1:3);
e2 = 1e3*err2(:,1:3);






