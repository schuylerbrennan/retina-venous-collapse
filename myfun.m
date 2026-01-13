function f = myfun(x,R,c,Q,L,mu,Pdown,IOP)

% wrong
% f = 2*(R+c*(Pdown + x/2 - IOP))*x^(1/4)-(128*mu*Q*L/pi)^(1/4);

% corrected
f = (((2*R + c*(Pdown + x/2 - IOP)))^4*x)-(128*mu*Q*L/pi);




