filename = "aceresults.txt";
fileID = fopen(filename,'r');
ace = fscanf(fileID, '%f');
fclose(fileID);
ace = reshape(ace,[3, 5])';

filename = "podresults.txt";
fileID = fopen(filename,'r');
pod = fscanf(fileID, '%f');
fclose(fileID);
pod = reshape(pod,[3, 5])';

filename = "pod23results.txt";
fileID = fopen(filename,'r');
pod23 = fscanf(fileID, '%f');
fclose(fileID);
pod23 = reshape(pod23,[3, 5])';

filename = "snapresults.txt";
fileID = fopen(filename,'r');
snap = fscanf(fileID, '%f');
fclose(fileID);
snap = reshape(snap,[3, 5])';

num=[5     5     5
    15    14    15
    31    30    31
    56    55    56
    92    91    91];
   
eerr = 1e3*[ace(:,2) snap(:,2) pod(:,2) pod23(:,2)];
ferr = 1e3*[ace(:,3) snap(:,3) pod(:,3) pod23(:,3)];

figure(1); clf; 
semilogy(num(:,1), eerr(:,1),'-o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on
semilogy(num(:,2), eerr(:,2), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,2), eerr(:,3), '--d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), eerr(:,4), '-.^', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'MarkerSize', 10); 
%semilogy(num(:,3), eerr(:,5), '-.v', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("Number of basis functions", 'FontSize', 18);
ylabel("MAE energy error (meV/atom)", 'FontSize', 18);
leg = legend({'ACE','SNAP-I','POD-I','POD-II'}, 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 95 0.025 1.2]);
set(gca, 'YTick', [0.03 0.06 0.12 0.24 0.48 0.96])
print -dpng GaN_train_energy_error.png

figure(2); clf; 
semilogy(num(:,1), ferr(:,1),'-o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on
semilogy(num(:,2), ferr(:,2), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,2), ferr(:,3), '--d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), ferr(:,4), '-.^', 'Color', [[0.4940 0.1840 0.5560]], 'LineWidth', 2, 'MarkerSize', 10); 
%semilogy(num(:,3), ferr(:,5), '-.v', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("Number of basis functions", 'FontSize', 18);
ylabel("MAE force error (meV/�)", 'FontSize', 18);
leg = legend({'ACE','SNAP-I','POD-I','POD-II'}, 'FontSize', 16, 'Location', 'NE');
% leg.ItemTokenSize = [50,10];
axis([5 95 6 80]);
set(gca, 'YTick', [7 10 14 20 28 40 56 80])
print -dpng GaN_train_force_error.png


% figure(2); clf; 
% plot(num(:,1), err(:,4), '-.o', 'LineWidth', 2, 'MarkerSize', 10); 
% hold on;
% plot(num(:,2), err(:,5), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
% plot(num(:,3), err(:,6), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
% set(gca,'FontSize',16); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% %title("Ta Eelement", 'FontSize', 18);
% xlabel("Number of basis functions", 'FontSize', 18);
% ylabel("MAE force error (meV/�)", 'FontSize', 18);
% legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 15, 'Location', 'NE');
% axis([5 150 0.04 0.38]);
% print -dpng InP_train_force_error.png

% filename = ['results/optpodlossfunction' num2str(5) '.txt'];
% fileID = fopen(filename,'r');
% frame = fscanf(fileID, '%f');
% fclose(fileID);
% sz = [41 41];
% a = reshape(frame,3,[])';
% X = reshape(a(:,1),sz);
% Y = reshape(a(:,2),sz);
% Z = reshape(a(:,3),sz);
% i1 = 1:41;
% i2 = 1:41;
% figure(3);clf; 
% pcolor(Y(i1,i2),X(i1,i2),Z(i1,i2));
% set(gca,'FontSize',16);
% set(gca,'LooseInset',get(gca,'TightInset'))
% [~,ind] = min(Z(:));
% u = [Y(ind) X(ind) Z(ind)];
% xlabel('Inner cut-off distance', 'FontSize',16);
% ylabel('Outer cut-off distance', 'FontSize',16);
% %hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% hold on; plot(u(1), u(2), 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFace', 'r');
% colorbar("FontSize",14);
% shading interp 
% print('-dpng', ['podlossfunction' num2str(5) '.png'])
% 
% filename = ['results/optacelossfunction' num2str(3) '.txt'];
% fileID = fopen(filename,'r');
% frame = fscanf(fileID, '%f');
% fclose(fileID);
% sz = [41 41];
% a = reshape(frame,3,[])';
% X = reshape(a(:,1),sz);
% Y = reshape(a(:,2),sz);
% Z = reshape(a(:,3),sz);
% i1 = 1:41;
% i2 = 1:41;
% figure(4);clf; 
% pcolor(Y(i1,i2),X(i1,i2),Z(i1,i2));
% set(gca,'FontSize',16);
% set(gca,'LooseInset',get(gca,'TightInset'))
% [~,ind] = min(Z(:));
% u = [Y(ind) X(ind) Z(ind)];
% xlabel('Inner cut-off distance', 'FontSize',16);
% ylabel('Outer cut-off distance', 'FontSize',16);
% %hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% hold on; plot(1.48, 4.5, 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFace', 'r');
% colorbar("FontSize",14);
% shading interp 
% print('-dpng', ['acelossfunction' num2str(3) '.png'])

% 
% Z = reshape(a(:,9),sz);
% figure(2);clf;
% mesh(X,Y,Z)
% [~,ind] = min(Z(:));
% u = [X(ind) Y(ind) Z(ind)];
% xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% zlabel('Force MAE', 'FontSize',12);
% hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% set(gca,'FontSize',14);
% print('-dpng', ['ForceMAE' num2str(n) '.png'])
% 
% filename = ['output' num2str(n) '.txt'];
% fileID = fopen(filename,'r');
% frame = fscanf(fileID, '%f');
% fclose(fileID);
% 
% sz = [37 55];
% a = reshape(frame,14,[])';
% X = reshape(a(:,1),sz);
% Y = reshape(a(:,2),sz);
% Z = reshape(a(:,5),sz);
% figure(3);clf;
% mesh(X,Y,Z)
% [~,ind] = min(Z(:));
% u = [X(ind) Y(ind) Z(ind)];
% xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% zlabel('Energy MAE', 'FontSize',12);
% hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% set(gca,'FontSize',14);
% print('-dpng', ['EnergyMAESurface' num2str(n) '.png'])
% 
% Z = reshape(a(:,5),sz);
% figure(4);clf;
% mesh(X,Y,Z)
% [~,ind] = min(Z(:));
% u = [X(ind) Y(ind) Z(ind)];
% xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% zlabel('Force MAE', 'FontSize',12);
% hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% set(gca,'FontSize',14);
% print('-dpng', ['ForceMAESurface' num2str(n) '.png'])
% 
% 
% 