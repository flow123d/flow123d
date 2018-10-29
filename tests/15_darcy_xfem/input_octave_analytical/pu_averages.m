function [press, vel, pr_avg, ur_avg] = pu_averages(xwa,ywa,rhoa,xw,yw, pr,ur)

% generate points on well edge
N = 1000;
Nw = length(xw);
x = zeros(N,1);
y = zeros(N,1);
phistep = 2*pi/N;
phi = 0;
for i = 1:N
  x(i) = rhoa * cos(phi) + xwa;
  y(i) = rhoa * sin(phi) + ywa;
  phi = phi + phistep;
endfor
L = 2*pi*rhoa;

%compute integrals
dl = 2*pi*rhoa / N;
press = zeros(1,Nw);
vel = zeros(1,Nw);
for w = 1:Nw
  if(xwa == xw(w) && ywa == yw(w))
    press(w) = log(rhoa) * L;
    vel(w) = 1/rhoa * L;
    continue;
  endif
  
  for j = 1:N
    rj = [x(j);y(j)] - [xw(w);yw(w)];
    ri = [x(j);y(j)] - [xwa;ywa];
    press(w) = press(w) + log(norm(rj)) * dl;
    vel(w) = vel(w) + (rj'*ri)/(rj'*rj*rhoa)*dl;
  endfor
endfor
press = press / L;
vel = vel / L;


% compute averages of regular functions
pr_avg = 0;
ur_avg = 0;
for i = 1:N
  pr_avg = pr_avg + pr(x(i),y(i)) * dl;
  n = [x(i)-xwa; y(i)-ywa] / rhoa;  % normal vector
  ur_avg = ur_avg + dot(ur(x(i),y(i)), n) * dl;
endfor
pr_avg = pr_avg / L;
ur_avg = ur_avg / L;

end