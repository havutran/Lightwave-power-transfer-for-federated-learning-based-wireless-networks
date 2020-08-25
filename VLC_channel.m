function [ h ] = VLC_channel( D1 )

theta = 60; % semi-angle at half power
ml=-log(2)/log(cosd(theta)); %Lambertian order of emission
Ar=85e-4; %detector physical area of a PD
Ts=1; %gain of an optical filter; ignore if no filter is used
index=1.5; %refractive index of a lens at a PD; ignore if no lens is used
FOV=70; %FOV of a receiver
G_Con=(index^2)/(sind(FOV).^2);
height=2; %the distance between source and receiver plane
%D1= 2.5; % distance vector from source 1
cosphi_A1 = height./D1;
receiver_angle=acosd(cosphi_A1);
if abs(receiver_angle)>FOV
    h =0;
else h = (ml+1)*Ar.*cosphi_A1.^(ml)./(2*pi.*D1.^2).*Ts.*G_Con.*cosd(receiver_angle);
end
%P_total = 6000;

end
