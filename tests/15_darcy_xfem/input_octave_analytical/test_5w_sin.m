clear all;
close all;
clc;

format long;

%data
K = 1e-3;
sigma = [10, 10, 10, 10, 10]';
rho = 0.03*[1, 1, 1, 1, 1]';
P = [-150, -30, 120, -50, 100]';

xw = [2.8, 4.9, 2.9, 7.3, 7.4]';
yw = [2.5, 5.4, 7.4, 7.8, 2.8]';

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

n = 200;
x = linspace(0,10,n);
y = linspace(0,10,n);
p = zeros(n);

for w = 1:length(a)
  for i = 1:n
    for j = 1:n
      r = [x(i)-xw(w); y(j)-yw(w)];
      p(i,j) = p(i,j) + a(w) * log(norm(r)) + U*sin(omg*x(j));
    endfor
  endfor
endfor

%ppi = @(x,y,i) a(i)*log(((x-xw(i))^2 + (y-yw(i))^2));
%p = @(x,y) ppi(x,y,1) + ppi(x,y,2) + ppi(x,y,3) + ppi(x,y,4) + ppi(x,y,5);
contour(x,y,p, 100)