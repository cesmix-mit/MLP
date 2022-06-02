a = 3.149;
b = 3.149*sqrt(3);
c = 3.480;
alpha = 90;
beta = 90;
gamma = 90;
atomtype = [1, 1, 2, 2, 2, 2];
atombasis =  [0   1.574500000000000   0.000000000000000   1.574500000000000   0.000000000000001   1.574500000000001;
              0   2.727113996517198   1.818075997678132   0.909037998839066   3.636151995356264   4.545189994195329;
              0                   0   1.740000000000000   1.740000000000000   1.740000000000000   1.740000000000000];
         
alpha = alpha*pi/180;
beta = beta*pi/180;
gamma = gamma*pi/180;
cosalpha = cos(alpha);
cosbeta = cos(beta);
cosgamma = cos(gamma);

lx = a;
xy = b*cosgamma;
xz = c*cosbeta;
ly = sqrt(b*b - xy*xy);
yz = (b*c*cosalpha - xy*xz)/ly;
lz = sqrt(c*c - xz*xz - yz*yz);

a1 = [lx; 0; 0];
a2 = [xy; ly; 0];
a3 = [xz; yz; lz];
fraccoords = inv([a1 a2 a3])*atombasis;

figure(3); clf;
plotboundingbox(a1, a2, a3);
plot3(atombasis(1,1),atombasis(2,1),atombasis(3,1),'or');    
plot3(atombasis(1,2),atombasis(2,2),atombasis(3,2),'or');    
plot3(atombasis(1,3),atombasis(2,3),atombasis(3,3),'ob');    
plot3(atombasis(1,4),atombasis(2,4),atombasis(3,4),'ob');    
plot3(atombasis(1,5),atombasis(2,5),atombasis(3,5),'ob');    
plot3(atombasis(1,6),atombasis(2,6),atombasis(3,6),'ob');    

m = 3;
n = 3;
p = 3;
xref = referenceorigins(m, n, p);
xlat = [a1 a2 a3]*xref;

figure(4); clf; hold on;
for i = 1:9
   plotboundingbox2(a1, a2, a3, xlat(:,i));
   plot3(atombasis(1,1)+xlat(1,i),atombasis(2,1)+xlat(2,i),atombasis(3,1)+xlat(3,i),'or');    
   plot3(atombasis(1,2)+xlat(1,i),atombasis(2,2)+xlat(2,i),atombasis(3,2)+xlat(3,i),'or');    
   plot3(atombasis(1,3)+xlat(1,i),atombasis(2,3)+xlat(2,i),atombasis(3,3)+xlat(3,i),'ob');    
   plot3(atombasis(1,4)+xlat(1,i),atombasis(2,4)+xlat(2,i),atombasis(3,4)+xlat(3,i),'ob');    
   plot3(atombasis(1,5)+xlat(1,i),atombasis(2,5)+xlat(2,i),atombasis(3,5)+xlat(3,i),'ob');    
   plot3(atombasis(1,6)+xlat(1,i),atombasis(2,6)+xlat(2,i),atombasis(3,6)+xlat(3,i),'ob');    
end


% i1 = (ti==1);
% plot3(xi(1,i1),xi(2,i1),xi(3,i1),'or');    
% i2 = (ti==2);
% plot3(xi(1,i2),xi(2,i2),xi(3,i2),'ob');    
% pause(0.1)
% if rem(i,100)==0
%     pause
% end






