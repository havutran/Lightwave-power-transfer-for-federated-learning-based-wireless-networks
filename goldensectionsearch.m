function [T] = goldensectionsearch (a,b,alpha,c,D,theta, Gamma, K, tau)
epsilon=0.000001;               % accuracy value
iter= 500;                       % maximum number of iterations
rho=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations
x1=a+(1-rho)*(b-a);             % computing x values
x2=a+rho*(b-a);

f_x1=fT(x1,alpha,c,D,theta, Gamma, K, tau);                     % computing values in x points
f_x2=fT(x2,alpha,c,D,theta, Gamma, K, tau);

while ((abs(b-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        b=x2;
        x2=x1;
        x1=a+(1-rho)*(b-a);
        
        f_x1=fT(x1,alpha,c,D,theta, Gamma, K, tau);
        f_x2=fT(x2,alpha,c,D,theta, Gamma, K, tau);
        
    else
        a=x1;
        x1=x2;
        x2=a+rho*(b-a);
        
        f_x1=fT(x1,alpha,c,D,theta, Gamma, K, tau);
        f_x2=fT(x2,alpha,c,D,theta, Gamma, K, tau);
    end
    
    k=k+1;
end

% chooses minimum point
if(f_x1<f_x2)
    T = x1;
else
    T = x2;

end