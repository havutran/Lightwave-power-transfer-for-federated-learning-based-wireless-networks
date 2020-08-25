function [P1] = bisection(P,H,epsilon,EH_VLC,tau)

fopt=0.75;
nu = 0.4;
I_d = 1e-9;
V= 26*10^(-3);

my_fun = @(x) (fopt*V*(log (nu*H*x/I_d + 1))*(nu*H*x) - (epsilon/tau - EH_VLC )) ;

a = 0;
b = P ;


while abs(a-b) > 1e-8
p = (a + b)/2;
if my_fun(p)<=0
    a=p;
else
    b=p;
end

end

P1 = p;

end