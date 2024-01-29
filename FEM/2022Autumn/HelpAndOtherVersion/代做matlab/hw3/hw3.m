clear;clc;close all;

q_ = -200;
F_ = -1000;
M_ = 2000;
L = 0.12;
d1 = 0.03;
d2 = 0.02;
Es = 200e9;
I1 = pi * d1^4 / 64;
I2 = pi * d2^4 / 64;
N = 51;
x = linspace(0, 2*L, N);
h = x(2) - x(1);

K = zeros(2*N, 2*N);
b = zeros(2*N);
localk = [12/h^3, 6/h^2, -12/h^3, 6/h^2;
   6/h^2, 4/h, -6/h^2, 2/h;
   -12/h^3, -6/h^2, 12/h^3, -6/h^2;
   6/h^2, 2/h, -6/h^2, 4/h];
localb = h * [1/2, h/12, 1/2, -h/12];
not_boundary = [3:2*N-2, 2*N];

for k = 1:N-1
    if k <= (N-1)/2
        K(2*k-1:2*k+2, 2*k-1:2*k+2) = K(2*k-1:2*k+2, 2*k-1:2*k+2) + Es * I1 * localk;
        b(2*k-1:2*k+2) = b(2*k-1:2*k+2) + q_ * localb;
    else
        K(2*k-1:2*k+2, 2*k-1:2*k+2) = K(2*k-1:2*k+2, 2*k-1:2*k+2) + Es * I2 * localk;
    end
end
b(N) = b(N) + F_;
b(2*N) = b(2*N) - M_;

bm = b(not_boundary);
Km = K(not_boundary, not_boundary);
solution = zeros(2*N);
solution(not_boundary) = Km \ bm.';
w = solution(1:2:2*N);
theta = solution(2:2:2*N);

plot(x, w)
xlabel('x')
ylabel('w')
title('w(x)-x')
exportgraphics(gca, 'w.png')
hold off;
plot(x, theta*180/pi)
xlabel('x')
ylabel('\theta')
title('\theta\circ(x)-x')
exportgraphics(gca, 'theta.png')
