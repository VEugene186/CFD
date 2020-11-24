clc;
clear all;
close all;

x = dlmread('x.csv', '\t');
y = dlmread('y.csv', '\t');
T = dlmread('T.csv', '\t');

f1 = figure(1);
colormap(hot(128));
contourf(x, y, T, 128, 'linecolor', 'none');
colorbar;
axis([-0.05, 1.05, -0.05, 1.05], 'square');