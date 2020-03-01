clc;
clear all;
close all;
x = dlmread('x.csv', '\t');
y = dlmread('y.csv', '\t');

plot(x, y, 'ro');