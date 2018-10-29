clear all;
clc;
close all;

format long ;

K = 1;
U = 80;
om = 1;
xw = 3.33;
yw = 3.33;
rho = 0.03;
sigma = 10;
pw = 100;

N = 10000;

x = zeros(N,1);
y = zeros(N,1);

% generate points on well edge
phistep = 2*pi/N;
phi = 0;
for i = 1:N
  x(i) = rho * cos(phi) + xw;
  y(i) = rho * sin(phi) + yw;
  phi = phi + phistep;
endfor
%phi
%2*pi
%plot(x,y);

%compute integral
dl = 2*pi*rho / N;
press = 0;
vel = 0;
for i = 1:N
  press = press + sin(om*x(i))*dl;
  vel = vel + K*om*cos(om*x(i))*(x(i)-xw)/rho*dl;
endfor
press = press / (2*pi*rho);
press_r_avg = U*press

vel = vel / (2*pi*rho);
vel_r_avg = U*vel

a = (press_r_avg - pw - vel_r_avg/sigma)/(K/sigma/rho - log(rho))

p_avg = a*log(rho) + press_r_avg
u_avg = K*a/rho + vel_r_avg

Ltest = sigma*(p_avg - pw)*2*pi*rho
Ptest = u_avg * 2*pi*rho

sum = 0;
%for i = 1:N
 % sum = sum + sigma(sin(om*x(i))*dl;
%endfor