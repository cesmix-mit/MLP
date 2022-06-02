function plotboundingbox(a, b, c)

w1 = [0; 0; 0];
w2 = a;
w3 = a+b;
w4 = b;
w5 = c;
w6 = a+c;
w7 = a+b+c;
w8 = b+c;

%figure(1); clf;
hold on;
plotline(w1, w2);
plotline(w2, w3);
plotline(w3, w4);
plotline(w4, w1);
plotline(w5, w6);
plotline(w6, w7);
plotline(w7, w8);
plotline(w8, w5);
plotline(w1, w5);
plotline(w2, w6);
plotline(w3, w7);
plotline(w4, w8);
axis equal
view(3);



