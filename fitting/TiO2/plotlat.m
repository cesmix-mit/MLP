filename = 'natom.bin';
fileID = fopen(filename,'r');
nat = fread(fileID,'double');
fclose(fileID);

filename = 'x.bin';
fileID = fopen(filename,'r');
x = fread(fileID,'double');
fclose(fileID);
x = reshape(x, 3, length(x)/3);

filename = 't.bin';
fileID = fopen(filename,'r');
t = fread(fileID,'double');
fclose(fileID);

filename = 'e.bin';
fileID = fopen(filename,'r');
e = fread(fileID,'double');
fclose(fileID);

filename = 'lattice.bin';
fileID = fopen(filename,'r');
lat = fread(fileID,'double');
fclose(fileID);
lat = reshape(lat, 9, length(lat)/9);

natom = [0; cumsum(nat)];

n = unique(nat);
ind = find(nat == 6);
nl = length(ind);
lv = zeros(nl,1);
le = zeros(nl,1);
lw = zeros(nl,1);
for i = 1:length(ind)
    ii = ind(i);
    a = lat(1:3,ii);
    b = lat(4:6,ii);
    c = lat(7:9,ii);
    le(i) = e(ii)/nat(ii);
    lv(i) = dot(a, cross(b,c))/nat(ii);
    lw(i) = dot(a, cross(b,c));
end

figure(1);clf; plot(lw/2,le-min(le),'o');
axis([20 55 0 1]);

lat6 = lat(:,ind);
a = [2 3 4 6 7 8];
m = max(abs(lat6(a,:)), [], 1);
i1 = find(m<1e-3);

ii = i1(411:500);
%figure(2);clf; plot(lv(ii),le(ii),'o');
figure(2);clf; plot(lw(ii)/2,le(ii)-min(le),'o');
axis([20 55 0 1]);
lw2 = lw(ii);
le2 = le(ii)-min(le);
ind2 = ind(ii);

lat6 = lat(:,ind);
a = [2 3 6 7 8];
m = max(abs(lat6(a,:)), [], 1);
i2 = find(m<1e-3);
i2 = setdiff(i2,i1);

ii = i2(1:end);
figure(3);clf; plot(lw(ii)/2,le(ii)-min(le),'o');
axis([20 55 0 1]);

lat6 = lat(:,ind);
a = [2 3 4 6];
m = max(abs(lat6(a,:)), [], 1);
i3 = find(m<1e-3);
i3 = setdiff(i3,union(i1,i2));

ii = i3(401:end);
figure(4);clf; plot(lw(ii)/2,le(ii)-min(le),'o');
axis([20 55 0 1]);
lw1 = lw(ii);
le1 = le(ii)-min(le);
ind1 = ind(ii);

i4 = setdiff(1:length(ind),union(union(i1,i2),i3));
ii = i4(1:end);
figure(5);clf; plot(lw(ii)/2,le(ii)-min(le),'o');
axis([20 55 0 1]);


ind = find(nat == 24);
nl = length(ind);
lv = zeros(nl,1);
le = zeros(nl,1);
lw = zeros(nl,1);
for i = 1:length(ind)
    ii = ind(i);
    a = lat(1:3,ii);
    b = lat(4:6,ii);
    c = lat(7:9,ii);
    le(i) = e(ii)/nat(ii);
    lv(i) = dot(a, cross(b,c))/nat(ii);
    lw(i) = dot(a, cross(b,c));
end

k = length(le);
ii = 2019:2109;
figure(6);clf; plot(lw(ii)/8,le(ii)-min(le),'o');
axis([20 55 0 1]);
lw3 = lw(ii);
le3 = le(ii)-min(le);
ind3 = ind(ii);

figure(7);clf; 
plot(lw1/2,le1,'o');
hold on;
plot(lw2/2,le2,'d');
plot(lw3/8,le3,'s');
axis([20 55 0 1]);


i = find(le1<=1);
x1 = lw1(i)/2;
y1 = le1(i);
[a,jj] = sort(x1);
x1 = x1(jj);
y1 = y1(jj);

i = find(le2<=1);
x2 = lw2(i)/2;
y2 = le2(i);
[a,jj] = sort(x2);
x2 = x2(jj);
y2 = y2(jj);

i = find(le3<=1);
x3 = lw3(i)/8;
y3 = le3(i);
[a,jj] = sort(x3);
x3 = x3(jj);
y3 = y3(jj);

figure(1);clf; 
plot(x1,y1,'o-k','LineWidth',2, 'MarkerSize', 8);
hold on;
plot(x2,y2,'s-', 'Color', [0.4660 0.6740 0.2880], 'LineWidth',2, 'MarkerSize', 8);
plot(x3,y3,'d-','Color', [0.5940 0.1840 0.5560], 'LineWidth',2, 'MarkerSize', 8);
axis([20 55 0 1]);
set(gca,'FontSize',18); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("Volume (Å^3/atom)", 'FontSize', 18);
ylabel("Energy (eV/atom)", 'FontSize', 18);
leg = legend({'Anatas', 'Rutile', 'Brookite'},'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
print -dpng TiO2_energy_curve.png


ind = find(nat == 24);
nl = length(ind);
lv = zeros(nl,1);
le = zeros(nl,1);
lw = zeros(nl,1);
for i = 1:length(ind)
    ii = ind(i);
    a = lat(1:3,ii);
    b = lat(4:6,ii);
    c = lat(7:9,ii);
    le(i) = e(ii)/nat(ii);
    lv(i) = dot(a, cross(b,c))/nat(ii);
    lw(i) = dot(a, cross(b,c));
end

y = le-min(le);
ij = find(y>1.15);
ii = ij(490:900);
x = lw(ii)/8;
y = le(ii)-min(le);

x1 = 40;
y1 = 1;
x2 = 70;
y2 = 2.7;
i = find( (x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0);
ii = ii(i);
x = lw(ii)/8;
y = le(ii)-min(le);

x1 = 30;
y1 = 1;
x2 = 25;
y2 = 1.8;
i = find( (x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0);
ii = ii(i);
x = lw(ii)/8;
y = le(ii)-min(le);

figure(2);clf; plot(x,y,'o');

x1 = 35;
y1 = 2.5;
x2 = 72;
y2 = 3.8;
x3 = 24;
y3 = 4.3;
i = find( ((x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0) & ((x-x1)/(x3-x1) - (y-y1)/(y3-y1) < 0));
s1 = x(i);
t1 = y(i);
j = setdiff(1:length(x),i);
x = x(j);
y = y(j);

figure(2);clf; plot(s1,t1,'o');
figure(3);clf; plot(x,y,'o');

x1 = 33;
y1 = 2.05;
x2 = 72;
y2 = 3.5;
x3 = 24;
y3 = 3.6;
i = find( ((x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0) & ((x-x1)/(x3-x1) - (y-y1)/(y3-y1) < 0));
s2 = x(i);
t2 = y(i);
j = setdiff(1:length(x),i);
figure(2);clf; plot(s2,t2,'o');

x = x(j);
y = y(j);
figure(3);clf; plot(x,y,'o');

x1 = 33;
y1 = 1.64;
x2 = 72;
y2 = 3.22;
x3 = 24;
y3 = 3.14;
i = find( ((x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0) & ((x-x1)/(x3-x1) - (y-y1)/(y3-y1) < 0));
s3 = x(i);
t3 = y(i);
j = setdiff(1:length(x),i);
figure(2);clf; plot(s3,t3,'o');

x = x(j);
y = y(j);
figure(3);clf; plot(x,y,'o');

x1 = 33.2;
y1 = 1.242;
x2 = 72;
y2 = 3.02;
x3 = 24;
y3 = 2.68;
i = find( ((x-x1)/(x2-x1) - (y-y1)/(y2-y1) < 0) & ((x-x1)/(x3-x1) - (y-y1)/(y3-y1) < 0));
s4 = x(i);
t4 = y(i);
j = setdiff(1:length(x),i);
figure(2);clf; plot(s4,t4,'o');

s5 = x(j);
t5 = y(j);
figure(3);clf; plot(s5,t5,'o');

[a,i3] = sort(s3); s3 = s3(i3); t3 = t3(i3);
[a,i3] = sort(s4); s4 = s4(i3); t4 = t4(i3);
[a,i3] = sort(s5); s5 = s5(i3); t5 = t5(i3);

figure(1);clf; 
plot(s5,t5,'o-k','LineWidth',2, 'MarkerSize', 8);
hold on;
plot(s3,t3,'s-', 'Color', [0.4660 0.6740 0.2880], 'LineWidth',2, 'MarkerSize', 8);
plot(s4,t4,'d-','Color', [0.5940 0.1840 0.5560], 'LineWidth',2, 'MarkerSize', 8);
axis([20 75 1 3.5]);
set(gca,'FontSize',18); 
set(gca,'LooseInset',max(get(gca,'TightInset'),0.005))
xlabel("Volume (Å^3/atom)", 'FontSize', 18);
ylabel("Energy (eV/atom)", 'FontSize', 18);
leg = legend({'Deformed Brookite 1', 'Deformed Brookite 2', 'Deformed Brookite 3'},'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [50,10];
print -dpng TiO2_energy_curve2.png










