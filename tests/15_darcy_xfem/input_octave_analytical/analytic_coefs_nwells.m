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
  A(i,i) = K/sigma(i)/rho(i);
  for j=1:n
    A(i,j) = A(i,j) - log(rho(j));
  end
  b(i) = -P(i);
end

A
b

alpha = A\b
