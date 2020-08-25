clear;
J = 3; %the number of users
M = 4; %the number of transmit antennas

%----------------------------------------------
%VLC channels
d1 = 2.3;
d2 = 2.2;
d3 = 2.1;
[ h01 ] = VLC_channel( d1 );
[ h02 ] = VLC_channel( d2 );
[ h03 ] = VLC_channel( d3 );

%VLC EH
[VLC_EH1] = VLC_EH (h01);
[VLC_EH2] = VLC_EH (h02);
[VLC_EH3] = VLC_EH (h03);

%----------------------------------------------
%IRL channels
[ h11 ] = IRL_channel( d1 );
[ h12 ] = IRL_channel( d2 );
[ h13 ] = IRL_channel( d3 );

%----------------------------------------------
%Communication and Computation
K = 5;
tau = 4;

theta = [30 35 40 45 50 55 60]*10^(3) ;%Kb -> b

S=0;S1=0;T1=0;T2=0;T3=0;Y1=0;Y2=0;Y3=0;

loop = 1000;

for l = 1:loop
%----------------------------------------------
%RF channels
Kdb = 8;
dist = 3;
[g1] = Pathloss_Rician_channels (dist+0.3,1,M, Kdb);
[g2] = Pathloss_Rician_channels (dist,1,M, Kdb);
[g3] = Pathloss_Rician_channels (dist-0.3,1,M, Kdb);

%Optimal beamforming
[w1,Gamma_1] = beamforming(g1,g2,g3);
[w2,Gamma_2] = beamforming(g2,g1,g3);
[w3,Gamma_3] = beamforming(g3,g2,g1);

%Solving SubOP1
alpha = 2*10^(-28);
c = 20;
D = 10*10^(6); %Mb -> bits
fmax = 1.5*10^(9); %GHz -> Hz
fmin = 0.3*10^(9); %GHz -> Hz
a = tau - c*D*K/fmin;
b = tau - c*D*K/fmax;

    for i = 1:length(theta)
            [T_trans1(i)] = goldensectionsearch (a,b,alpha,c,D,theta(i), Gamma_1, K, tau);
            [T_trans2(i)] = goldensectionsearch (a,b,alpha,c,D,theta(i), Gamma_2, K, tau);
            [T_trans3(i)] = goldensectionsearch (a,b,alpha,c,D,theta(i), Gamma_3, K, tau);

            %Optimal Pu
            B = 10^(6);
            Pu1(i) = (2^(theta(i)/(T_trans1(i)*B)) - 1)/Gamma_1;

            Pu2(i) = (2^(theta(i)/(T_trans2(i)*B)) - 1)/Gamma_2;

            Pu3(i) = (2^(theta(i)/(T_trans3(i)*B)) - 1)/Gamma_3;

            %Optimal P1j
            P = 120;

            epsilon1(i) = fT(T_trans1(i),alpha,c,D,theta(i), Gamma_1, K, tau);
            epsilon2(i) = fT(T_trans2(i),alpha,c,D,theta(i), Gamma_2, K, tau);
            epsilon3(i) = fT(T_trans3(i),alpha,c,D,theta(i), Gamma_3, K, tau);

        for j = 1:length(h11)
            [P11(j,i)] = bisection(P,h11(j),epsilon1(i),VLC_EH1,tau);
            [P12(j,i)] = bisection(P,h12(j),epsilon2(i),VLC_EH2,tau);
            [P13(j,i)] = bisection(P,h13(j),epsilon3(i),VLC_EH3,tau);
        end
    end
    S = S + P11 + P12 + P13;
    S1 = S1 + P11 + P12;
    T1 = T1 + T_trans1(3);
    T2 = T2 + T_trans2(3);
    T3 = T3 + T_trans3(3);
    
    Y1 = Y1 + T_trans1(6);
    Y2 = Y2 + T_trans2(6);
    Y3 = Y3 + T_trans3(6);
end

P1 = S/loop;
P2 = S1/loop;
T1 = T1/loop;
T2 = T2/loop;
T3 = T3/loop;
Y1 = Y1/loop;
Y2 = Y2/loop;
Y3 = Y3/loop;

fig1 = figure(1);
plot(theta*10^(-3),P1(1,:), '--r', theta*10^(-3),P1(2,:), ':ob', theta*10^(-3),P2(1,:),'-k',theta*10^(-3),P2(2,:),'-.*g','linewidth',2)
grid on
ylabel('Total Transmit Power of IRL (W)')
xlabel('Uplink Rate Achieved at Each User (Kbps)')
legend('3 users - \phi_{1,1/2} = 6^{\circ}', '3 users - \phi_{1,1/2} = 8^{\circ}', '2 users- \phi_{1,1/2} = 6^{\circ}', '2 users - \phi_{1,1/2} = 8^{\circ}')
saveas(fig1,'Power_r1', 'epsc')


figure(2);
hold off
subplot(1,2,1)
bar([1 2 3], tau*[1 1 1],'b' )
hold on
bar([1 2 3], [T1 T2 T3],'g' )
grid on
title('\theta = 40 kbps')
xlabel('User Index (n)')
ylabel('Time')
legend('Computation Time','Transmission Time')
subplot(1,2,2)
bar([1 2 3], tau*[1 1 1],'b' )
hold on
bar([1 2 3], [Y1 Y2 Y3],'g' )
grid on
title('\theta = 55 kbps')
xlabel('User Index (n)')
ylabel('Time')
legend('Computation Time','Transmission Time')

c = categorical({'1','2','3','1','2','3'});






