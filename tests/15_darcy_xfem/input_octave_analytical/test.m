clear all;
close all;
clc;

rho = 0.03
sigma = 10

run "output/matrix.m"

spy(Mat_0x84000000_0);

s = Vec_0x84000000_2;

s(57)
2*pi*rho*sigma*(s(60) - s(113))



