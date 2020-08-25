function [Hpl] = Pathloss_Rician_channels (dist,n_u,n_t, Kdb)
pl = 2.6;

K = 10.^(Kdb/10); %Rician factor

%Rician channels
mu = sqrt( K/((K+1))); 

s = sqrt( 1/(2*(K+1)));

Hw = mu + s*(randn(n_t,n_u) + 1j*randn(n_t,n_u));

Hpl = [sqrt(1./(dist(1)^pl))*ones(n_t,1) ].*Hw;
end