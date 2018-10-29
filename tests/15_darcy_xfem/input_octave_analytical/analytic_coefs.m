function [alpha] = analytic_coefs(K,sigma,rho,P,xw,yw,source_p, source_u);

n = length(sigma);
A = zeros(n,n);
b = zeros(n,1);


for i=1:n
  [logr_avg, r_avg,  pr, unr] = pu_averages(xw(i), yw(i), rho(i), xw, yw, source_p, source_u);
  A(i,:) = A(i,:) + K*r_avg - sigma(i)*logr_avg;
  b(i) = sigma(i)*(pr-P(i)) - unr;
end

%A
%b
alpha = A\b;
