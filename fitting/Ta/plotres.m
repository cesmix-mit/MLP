filename = "numdesc.txt";
fileID = fopen(filename,'r');
num = fscanf(fileID, '%f');
fclose(fileID);

num = reshape(num,[3, 6])';

filename = "trainerror.txt";
fileID = fopen(filename,'r');
err = fscanf(fileID, '%f');
fclose(fileID);
err = reshape(err,[6, 6])';

err(:,1:3) = 1e3*err(:,1:3);
figure(1); clf; 
semilogy(num(:,1), err(:,1),'-.o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on
semilogy(num(:,2), err(:,2), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
semilogy(num(:,3), err(:,3), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
%title("Ta Eelement", 'FontSize', 18);
xlabel("Number of descriptors", 'FontSize', 18);
ylabel("MAE energy error (meV/atom)", 'FontSize', 18);
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 150 0.9 500]);
print -dpng Ta_train_energy_error.png

figure(2); clf; 
plot(num(:,1), 1e3*err(:,4), '-.o', 'LineWidth', 2, 'MarkerSize', 10); 
hold on;
plot(num(:,2), 1e3*err(:,5), '--s', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10); 
plot(num(:,3), 1e3*err(:,6), '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10); 
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
%title("Ta Eelement", 'FontSize', 18);
xlabel("Number of descriptors", 'FontSize', 18);
ylabel("MAE force error (meV/Å)", 'FontSize', 18);
leg = legend({'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 150 0.04*1e3 0.38*1e3]);
print -dpng Ta_train_force_error.png

fcc = [ 6.750000000000000  10.085400000000000  10.086487297023403  10.085583112705418  10.085482526872017
   7.447750000000000   4.904480000000000   4.894688830945612   4.903789897394369   4.904014007420820
   8.192000000000000   0.865768000000000   0.887451877024437   0.867646692138813   0.866877659114349
   8.984249999999999  -2.291170000000000  -2.293988230748911  -2.294240027849480  -2.292611556720869
   9.826000000000001  -4.769000000000000  -4.791783661028304  -4.765964880309853  -4.768245695992498
  10.718800000000000  -6.715680000000000  -6.725396439354384  -6.716520315045134  -6.714960864582323
  11.664000000000000  -8.229939999999999  -8.213836096792230  -8.232372755485375  -8.231720581007211
  12.663300000000000  -9.386860000000000  -9.367404675804407  -9.384435905742899  -9.385940454801169
  13.718000000000000 -10.247999999999999 -10.239029405210818 -10.247193104738546 -10.245570118932646
  14.829700000000001 -10.863500000000000 -10.867321722133646 -10.867922242065802 -10.868500926354210
  16.000000000000000 -11.275100000000000 -11.282978241680622 -11.273237456792714 -11.274158763829591
  17.230200000000000 -11.513100000000000 -11.517900262961001 -11.512577999478088 -11.511722203398351
  18.521999999999998 -11.606700000000000 -11.609594274296962 -11.609368925842466 -11.609352007113779
  19.876700000000000 -11.580900000000000 -11.584024335390337 -11.581571033504611 -11.584062258198720
  21.295999999999999 -11.457200000000000 -11.461260682536519 -11.455570559439662 -11.456327802291794
  22.781199999999998 -11.254000000000000 -11.258686264176532 -11.253350369854067 -11.251183209680692
  24.334000000000000 -10.987500000000001 -10.991513231469977 -10.989071245047226 -10.987213446271836
  25.955800000000000 -10.671500000000000 -10.673557935071122 -10.672328076671473 -10.672087236822406
  27.648000000000000 -10.318199999999999 -10.317901868939320 -10.319172184080934 -10.318246672342680
  29.412299999999998  -9.938330000000001  -9.937291224621328  -9.940029004004943  -9.938747021871215
  31.250000000000000  -9.540070000000000  -9.544266172100086  -9.540458156088839  -9.541069060999273
  33.162700000000001  -9.131489999999999  -9.142868995489874  -9.130672568770008  -9.131918517085881
  35.152000000000001  -8.718400000000001  -8.733279621961142  -8.719796474945605  -8.717945000135421
  37.219200000000001  -8.306160000000000  -8.321325288333359  -8.309609119213727  -8.305331676498218
  39.366000000000000  -7.898990000000000  -7.912010421672214  -7.901127901802274  -7.898907433771459
  41.593800000000002  -7.500300000000000  -7.509668116283130  -7.500707246151273  -7.501595737079103
  43.904000000000003  -7.112790000000000  -7.118071602442957  -7.114377770204545  -7.114764981659635
  46.298299999999998  -6.738600000000000  -6.740514990167051  -6.741853139793529  -6.739457146989874
  48.777999999999999  -6.379970000000000  -6.379872288343090  -6.380773539574613  -6.377945488936481
  51.344799999999999  -6.037590000000000  -6.038642290790485  -6.034991006885264  -6.034919036807511
  54.000000000000000  -5.712600000000000  -5.718985102616624  -5.717925316416610  -5.717845484156562];

figure(3); clf; 
plot(fcc(:,1),fcc(:,2),'o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(fcc(:,1),fcc(:,3), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10);
plot(fcc(:,1),fcc(:,4), '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10);
plot(fcc(:,1),fcc(:,5), '-.', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'MarkerSize', 10);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
%title("Ta Eelement", 'FontSize', 18);
xlabel("Volume (Å^3/atom)", 'FontSize', 17);
ylabel("Energy (eV/atom)", 'FontSize', 17);
leg = legend({'DFT', 'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 55 -12 12]);
print -dpng Ta_energy_curve_fcc.png


bcc=[5.324000000000000  28.131399999999999  28.131018172586636  28.131411949057565  28.131375239493796
   6.083500000000000  16.555000000000000  16.558279008148968  16.554988777766432  16.555028420457766
   6.912000000000000   8.157069999999999   8.146623089929781   8.157104415713729   8.157183353580731
   7.812500000000000   2.048440000000000   2.063354531162445   2.048424460830566   2.048108820748170
   8.788000000000000  -2.390660000000000  -2.399339352784744  -2.390524712253148  -2.390287804688914
   9.841500000000000  -5.598600000000000  -5.598773880802095  -5.599801291738529  -5.598502614519888
  10.976000000000001  -7.890260000000000  -7.884372097767278  -7.887053639509520  -7.890783907940776
  12.194500000000000  -9.501099999999999  -9.504473744303763  -9.506609402957857  -9.499810889167408
  13.500000000000000 -10.604500000000000 -10.611648585753573 -10.597802109269239 -10.604797054772041
  14.895500000000000 -11.310300000000000 -11.314905579255841 -11.313202315052894 -11.309919541019868
  16.384000000000000 -11.702700000000000 -11.705112517640742 -11.706613676279057 -11.700583956777784
  17.968499999999999 -11.848400000000000 -11.843395487574226 -11.842141223360057 -11.841761496521544
  19.652000000000001 -11.800700000000001 -11.786910340868493 -11.796159593970849 -11.795031496146089
  21.437500000000000 -11.603700000000000 -11.586820513101364 -11.603787185649614 -11.603690365755529
  23.327999999999999 -11.292600000000000 -11.284951024232004 -11.292889095482449 -11.294896418124054
  25.326499999999999 -10.897800000000000 -10.896875712460552 -10.900205898579792 -10.897700641368754
  27.436000000000000 -10.443000000000000 -10.446547047377512 -10.444021409538216 -10.441543092680961
  29.659500000000001  -9.949040000000000  -9.952992476879686  -9.946360965663112  -9.948689479739922
  32.000000000000000  -9.432470000000000  -9.430895270799660  -9.433267304886597  -9.433188529988092
  34.460500000000003  -8.906540000000000  -8.892521233628734  -8.912703515097956  -8.906643507304512
  37.043999999999997  -8.381810000000000  -8.348921462158101  -8.382861151892733  -8.381836760760452];

figure(4); clf; 
plot(bcc(:,1),bcc(:,2),'o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
plot(bcc(:,1),bcc(:,3), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,4), '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'MarkerSize', 10);
plot(bcc(:,1),bcc(:,5), '-.', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2, 'MarkerSize', 10);
set(gca,'FontSize',16); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("Volume (Å^3/atom)", 'FontSize', 17);
ylabel("Energy (eV/atom)", 'FontSize', 17);
leg = legend({'DFT', 'ACE','SNAP','POD'},'interpreter', 'latex', 'FontSize', 16, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
axis([5 38 -12.5 30]);
print -dpng Ta_energy_curve_bcc.png


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
% % 
% % Z = reshape(a(:,9),sz);
% % figure(2);clf;
% % mesh(X,Y,Z)
% % [~,ind] = min(Z(:));
% % u = [X(ind) Y(ind) Z(ind)];
% % xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% % ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% % zlabel('Force MAE', 'FontSize',12);
% % hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% % set(gca,'FontSize',14);
% % print('-dpng', ['ForceMAE' num2str(n) '.png'])
% % 
% % filename = ['output' num2str(n) '.txt'];
% % fileID = fopen(filename,'r');
% % frame = fscanf(fileID, '%f');
% % fclose(fileID);
% % 
% % sz = [37 55];
% % a = reshape(frame,14,[])';
% % X = reshape(a(:,1),sz);
% % Y = reshape(a(:,2),sz);
% % Z = reshape(a(:,5),sz);
% % figure(3);clf;
% % mesh(X,Y,Z)
% % [~,ind] = min(Z(:));
% % u = [X(ind) Y(ind) Z(ind)];
% % xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% % ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% % zlabel('Energy MAE', 'FontSize',12);
% % hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% % set(gca,'FontSize',14);
% % print('-dpng', ['EnergyMAESurface' num2str(n) '.png'])
% % 
% % Z = reshape(a(:,5),sz);
% % figure(4);clf;
% % mesh(X,Y,Z)
% % [~,ind] = min(Z(:));
% % u = [X(ind) Y(ind) Z(ind)];
% % xlabel('Cut-off radius', 'FontSize',12, 'Rotation', 18);
% % ylabel('Scaling parameter', 'FontSize',12, 'Rotation', -26);
% % zlabel('Force MAE', 'FontSize',12);
% % hold on; plot3(u(1), u(2), u(3), '*r', 'MarkerSize', 10, 'LineWidth', 2);
% % set(gca,'FontSize',14);
% % print('-dpng', ['ForceMAESurface' num2str(n) '.png'])
% % 
% % 
% % 
