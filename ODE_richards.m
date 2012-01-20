clear all;
clear global;

global n=2;
global m=1-1/n;
global h_d =1;
global dt =0.1;
global h0 =-10;
tend=100;

function y = th(h)
  global n;
  global m;	
  if (h<0)
  	y=(1+(-h)**n)**(-m);
   else
  	y=1;
   endif
endfunction

function y = cap(h)
  global n;
  global m;	
  if (h<0)
  	y=(1 + (-h)**n)**(-2. + 1./n) * (-h)**(n-1) *(n-1);
   else
  	y=1;
   endif
endfunction


function y = con(h)
  global n;
  global m;	
  if (h<0)
  	y=th(h)**(0.5) * ( 1 - ( 1 - th(h)**(1/m) )**(m))**2;
   else
  	y=1;
   endif
endfunction

function y = f(h)
global h_d;
global dt;
global h_last;
   c=con((h+h_last)/2);
   y=0.01*(th(h)-th(h_last))/dt+(0.5*c*(h-h_d)+0.5*c*(h_last-h_d) );
endfunction

global h_last;
i=0;
h_last=h0;      
nn=round(tend/dt);
x=(1:nn);
y=x;


for t=0:dt:tend
  x(i+1)=t;
  y(i+1) = fzero(@f, [-10,100]);  %,optimset("TolX",1e-10));
  i=i+1;
  h_last = y(i);
endfor	



% y=arrayfun(@con,h);
plot(x,y)
