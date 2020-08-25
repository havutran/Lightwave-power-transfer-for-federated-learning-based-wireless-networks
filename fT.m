function y = fT(T,alpha,c,D,theta, Gamma, K, tau)
B = 10^(6);
y = ((alpha/2)*(c*D)^3*K^2)/(tau-T)^2 + T*( 2^(theta/(T*B)) -1)/Gamma;
end