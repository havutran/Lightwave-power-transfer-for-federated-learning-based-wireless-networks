function [EH] = IVLC_EH (D1)

f=0.75;
nu = 0.4;
I_D = 1e-9;
V_0= 26*10^(-3);
P_vlc = 30*7*3*3; %total power of each LED array
P_sel = P_vlc; %equal-element equal-color power 


I_H = 12*10^(-3);
I_L = 2*10^(-3);

h = IVLC_channel(D1);

V = V_0*log ((nu*h*(P_vlc*I_H))/I_D + 1);

EH = f*( (nu*h*(P_vlc*I_H ))  )*V;

end