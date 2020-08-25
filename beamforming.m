function [w, Gamma] = beamforming(g1,g2,g3)
M = length(g1);
sigma2 = 10^(-10);
[A,B] = eig(g1*g1'*(g2*g2' + g3*g3' + sigma2*eye(M))^(-1));
[X,I] = max(diag(B));
w =  A(:,I);
Gamma = abs((w'*g1*g1'*w)/(w'*(g2*g2' + g3*g3' + sigma2*eye(M))*w ));
end