clear;clc;close all;

NE = textread('mesh.txt', '%d', 2);
N = NE(1);
E = NE(2);
data = readmatrix('mesh.txt');
position = data(1:N, 1:2);
connect = round(data(N+1:N+E, 2:4)) + 1;
b = 10;
alpha = 30*pi/180;

NI = [];
toler = 1e-5;
for i = 1: N
    x = position(i, 1);
    y = position(i, 2);
    a1 = abs(x - b);
    a2 = abs(y - tan(alpha) * x);
    a3 = abs(y + tan(alpha) * x);
    if a1 > toler && a2 > toler && a3 > toler
        NI = [NI, i];
    end
end

K = zeros(N, N);
F = zeros(1, N);
u = zeros(1, N);
for k = 1: E
    node1 = connect(k, 1);
    node2 = connect(k, 2);
    node3 = connect(k, 3);
    x1 = position(node1, 1);
    y1 = position(node1, 2);
    x2 = position(node2, 1);
    y2 = position(node2, 2);
    x3 = position(node3, 1);
    y3 = position(node3, 2);
    A2 = 2*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
    k11 = (x2 - x3)^2 + (y2 - y3)^2;
    k12 = (x1 - x3) * (-x2 + x3) + (y1 - y3) * (-y2 + y3);
    k13 = (x1 - x2) * (x2 - x3) + (y1 - y2) * (y2 - y3);
    k22 = (x1 - x3)^2 + (y1 - y3)^2;
    k23 = -(x1^2 + x2 * x3 - x1 * (x2 + x3) + (y1 - y2) * (y1 - y3));
    k33 = (x1 - x2)^2 + (y1 - y2)^2;
    F1 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3)) / 6;
    F2 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1+y3)) / 6;
    F3 = (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1+y3)) / 6;
    K(node1, node1) = K(node1, node1) + k11 / A2;
    K(node1, node2) = K(node1, node2) + k12 / A2;
    K(node1, node3) = K(node1, node3) + k13 / A2;
    K(node2, node1) = K(node2, node1) + k12 / A2;
    K(node2, node2) = K(node2, node2) + k22 / A2;
    K(node2, node3) = K(node2, node3) + k23 / A2;
    K(node3, node1) = K(node3, node1) + k13 / A2;
    K(node3, node2) = K(node3, node2) + k23 / A2;
    K(node3, node3) = K(node3, node3) + k33 / A2;
    F(node1) = F(node1) + F1;
    F(node2) = F(node2) + F2;
    F(node3) = F(node3) + F3;
end

u(NI) = K(NI, NI) \ F(NI).';
trisurf(connect, position(:, 1), position(:, 2), u)

view([-62.63 75.27])
exportgraphics(gca, 'res1.png')
view([-0.44 90.00])
exportgraphics(gca, 'res2.png')
view([55.41 66.95])
exportgraphics(gca, 'res3.png')
