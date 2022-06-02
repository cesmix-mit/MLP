filename = "results/ac3.txt";
fileID = fopen(filename,'r');
tm = fscanf(fileID, '%f');
fclose(fileID);
n = length(tm)/2;
ac3 = reshape(tm,[2, n])';

filename = "results/energyg3.txt";
fileID = fopen(filename,'r');
e3 = fscanf(fileID, '%f');
fclose(fileID);

sz = [25 25];
X = reshape(ac3(:,1),sz);
Y = reshape(ac3(:,2),sz);
Z = reshape(e3,sz);
i1 = 1:25;
i2 = 1:25;
[Zmin,ind] = min(Z(:));
Z = Z - Zmin;
u = [X(ind) Y(ind) Z(ind)];

figure(1);clf; 
pcolor(X(i1,i2),Y(i1,i2),Z(i1,i2));
set(gca,'FontSize',18);
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel('$a$', 'interpreter','latex', 'FontSize',26);
ylabel('$c$', 'interpreter','latex', 'FontSize',26);
hold on; plot(u(1), u(2), 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFace', 'r');
colorbar("FontSize",18);
shading interp 
print('-dpng', 'DFTenergy3.png')


filename = "results/fitsnapestenergies4.txt";
fileID = fopen(filename,'r');
tm = fscanf(fileID, '%f');
fclose(fileID);
n = length(tm)/2;
e = reshape(tm,[2, n])';
e = e(:,1) - Zmin;
e = reshape(abs(e(:)-Z(:)),sz);

figure(2);clf; 
pcolor(X(i1,i2),Y(i1,i2),e(i1,i2));
set(gca,'FontSize',18);
xlabel('$a$', 'interpreter','latex', 'FontSize',26);
ylabel('$c$', 'interpreter','latex', 'FontSize',26);
colorbar("FontSize",18);
shading interp 
axis tight;
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.04))
print('-dpng', 'SNAPenergy3.png')

filename = "results/fitsnapestenergies4.txt";
fileID = fopen(filename,'r');
tm = fscanf(fileID, '%f');
fclose(fileID);
n = length(tm)/2;
e = reshape(tm,[2, n])';
e = e(:,1) - Zmin;
e = reshape(e,sz);

filename = "results/fitpod23testenergies4.txt";
fileID = fopen(filename,'r');
tm = fscanf(fileID, '%f');
fclose(fileID);
n = length(tm)/2;
e23 = reshape(tm,[2, n])';
e23 = e23(:,1) - Zmin;
e23 = reshape(e23,sz);

filename = "results/fitaceestenergies4.txt";
fileID = fopen(filename,'r');
tm = fscanf(fileID, '%f');
fclose(fileID);
n = length(tm)/2;
eace = reshape(tm,[2, n])';
eace = eace(:,1) - Zmin;
eace = reshape(eace,sz);

j = 1;
bcc = [X(:,j) Z(:,j) e(:,j) eace(:,j) e23(:,j)];
figure(4); clf; 
plot(bcc(:,1),bcc(:,2),'o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(bcc(:,1),bcc(:,3), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,4), '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,5), '-.k', 'LineWidth', 2, 'MarkerSize', 10);
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel('$a$', 'interpreter','latex', 'FontSize',24);
ylabel("Energy (eV)", 'FontSize', 20);
leg = legend({'DFT', 'SNAP','ACE','POD-II'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([2.6 3.7 0 7]);
print -dpng HfB2_energy_curve.png

j = 10;
bcc = [X(:,j) Z(:,j) e(:,j) eace(:,j) e23(:,j)];
figure(4); clf; 
plot(bcc(:,1),bcc(:,2),'o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(bcc(:,1),bcc(:,3), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,4), '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,5), '-.k', 'LineWidth', 2, 'MarkerSize', 10);
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel('$a$', 'interpreter','latex', 'FontSize',24);
ylabel("Energy (eV)", 'FontSize', 20);
leg = legend({'DFT', 'SNAP','ACE','POD-II'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([2.6 3.7 0 4.2]);
print -dpng HfB2_energy_curve2.png

j = 25;
bcc = [X(:,j) Z(:,j) e(:,j) eace(:,j) e23(:,j)];
figure(4); clf; 
plot(bcc(:,1),bcc(:,2),'o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(bcc(:,1),bcc(:,3), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,4), '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,5), '-.k', 'LineWidth', 2, 'MarkerSize', 10);
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel('$a$', 'interpreter','latex', 'FontSize',24);
ylabel("Energy (eV)", 'FontSize', 20);
leg = legend({'DFT', 'SNAP','ACE','POD-II'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'SE');
leg.ItemTokenSize = [50,10];
axis([2.6 3.7 0 4]);
print -dpng HfB2_energy_curve3.png


% filename = "results/ac6.txt";
% fileID = fopen(filename,'r');
% tm = fscanf(fileID, '%f');
% fclose(fileID);
% n = length(tm)/2;
% ac3 = reshape(tm,[2, n])';
% 
% filename = "results/energyg6.txt";
% fileID = fopen(filename,'r');
% e3 = fscanf(fileID, '%f');
% fclose(fileID);
% 
% sz = [25 25];
% X = reshape(ac3(:,1),sz);
% Y = reshape(ac3(:,2),sz);
% Z = reshape(e3,sz);
% i1 = 1:25;
% i2 = 1:25;
% [Zmin,ind] = min(Z(:));
% Z = Z - Zmin;
% u = [X(ind) Y(ind) Z(ind)];
% 
% figure(2);clf; 
% pcolor(X(i1,i2),Y(i1,i2),Z(i1,i2));
% set(gca,'FontSize',18);
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel('$a$', 'interpreter','latex', 'FontSize',26);
% ylabel('$c$', 'interpreter','latex', 'FontSize',26);
% hold on; plot(u(2), u(1), 'or', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFace', 'r');
% colorbar("FontSize",18);
% shading interp 
% print('-dpng', 'DFTenergy6.png')
% 

