clc;
clear all;
close all;

nXMin = 5;
nYMin = 8;
nXMax = 13;
nYMax = 12;

dx = nXMax - nXMin;
dy = nYMax - nYMin;
N = 2 * (dx + dy);

cornerRightTop = 1 + dy / 2;
cornerLeftTop = cornerRightTop + dx;
cornerLeftBottom = cornerLeftTop + dy;
cornerRightBottom = cornerLeftBottom + dx;

t = linspace(0, 2 * pi, N + 1);
t = t(1 : N);

t = [t(cornerLeftBottom + 1: end) , t(1 : cornerLeftBottom)]

x = 0.5 + 0.1 * cos(t);
y = 0.5 + 0.1 * sin(t);
plot(x, y, 'ro');
axis([0 1 0 1], 'square');
grid on;
hold on;
plot(x(end), y(end), 'kx');
plot(x(end - dy), y(end - dy), 'k+');
%plot(x(cornerLeftBottom), y(cornerLeftBottom), 'k>');
%plot(x(cornerRightBottom), y(cornerRightBottom), 'k<');

dlmwrite('profile.dat', [ fliplr(x).' , fliplr(y).' ], '\t');