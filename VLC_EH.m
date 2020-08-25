function [VLC_EH] = VLC_EH (h)

f=0.75;
nu = 0.4;
I_d = 1e-9;
V_0= 26*10^(-3);

P_vlc = 28; 

V = V_0*log ((nu*h*P_vlc)/I_d + 1);

VLC_EH = f*( (nu*h*P_vlc)  )*V;

end