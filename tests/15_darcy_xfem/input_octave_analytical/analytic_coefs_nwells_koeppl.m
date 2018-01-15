clear all;
close all;
clc;

format long;

%data
K = 1.0;
sigma = [10, 10]';
rho = [0.03, 0.03]';
P = [50, 50]';

x = [4.42, 2.08]';
y = [2.08, 4.42]';

n = length(sigma);
A = zeros(n,n);
b = zeros(n,1);


for i=1:n
  A(i,i) = 1.0 + sigma(i);
  for j=1:n
    if i ~= j
      r = (x(i)-x(j))*(x(i)-x(j)) + (y(i)-y(j))*(y(i)-y(j))
      r = sqrt(r);
      A(i,j) = A(i,j) + sigma(i)*(1-rho(j)*log(r/rho(j)));
    endif
  end
  b(i) = sigma(i) * P(i);
end

A
b

alpha = A\b
