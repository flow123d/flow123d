clear all;
close all;
clc;

format long;

%data
K = 1.0;
sigma = [10, 10]';
rho = [0.03, 0.03]';
P = [100, 80]';

xw = [4.42, 2.08]';
yw = [2.08, 4.42]';

n = length(sigma);
A = zeros(n,n);
b = zeros(n,1);

%r_avg = zeros(n,n);
%logr_avg = zeros(n,n);

U = 80;
omg = 1;

source_p = @(x,y) U*sin(omg*x);
source_u = @(x,y) U*omg*cos(omg*x)*[1; 0];

a = analytic_coefs(K,sigma,rho,P,xw,yw,source_p, source_u)
