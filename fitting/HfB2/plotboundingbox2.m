function plotboundingbox2(a, b, c, w)

w1 = w;
w2 = w+a;
w3 = w+a+b;
w4 = w+b;
w5 = w+c;
w6 = w+a+c;
w7 = w+a+b+c;
w8 = w+b+c;

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



