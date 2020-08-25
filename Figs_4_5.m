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
K = [1 2 3 4 5];
tau = [3 4];

theta = 40*10^(3) ;%Kb -> b

S1=0;S2=0;S3=0;T1=0;T2=0;T3=0;

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


%plot K vs. Tcomp with various tau
    for i = 1:length(K)
        for j = 1:length(tau)
            a = tau(j) - c*D*K(i)/fmin;
            b = tau(j) - c*D*K(i)/fmax;
            [T_trans1(j,i)] = goldensectionsearch (a,b,alpha,c,D,theta, Gamma_1, K(i), tau(j));
            [T_trans2(j,i)] = goldensectionsearch (a,b,alpha,c,D,theta, Gamma_2, K(i), tau(j));
            [T_trans3(j,i)] = goldensectionsearch (a,b,alpha,c,D,theta, Gamma_3, K(i), tau(j));

            %Optimal Pu
            B = 10^(6);
            Pu1(j,i) = (2^(theta/(T_trans1(j,i)*B)) - 1)/Gamma_1;

            Pu2(j,i) = (2^(theta/(T_trans2(j,i)*B)) - 1)/Gamma_2;

            Pu3(j,i) = (2^(theta/(T_trans3(j,i)*B)) - 1)/Gamma_3;

            %Optimal P1j
            P = 120;

            epsilon1(j,i) = fT(T_trans1(j,i),alpha,c,D,theta, Gamma_1, K(i), tau(j));
            epsilon2(j,i) = fT(T_trans2(j,i),alpha,c,D,theta, Gamma_2, K(i), tau(j));
            epsilon3(j,i) = fT(T_trans3(j,i),alpha,c,D,theta, Gamma_3, K(i), tau(j));
        

            [P11(j,i)] = bisection(P,h11(2),epsilon1(j,i),VLC_EH1, tau(j));
            [P12(j,i)] = bisection(P,h12(2),epsilon2(j,i),VLC_EH2, tau(j));
            [P13(j,i)] = bisection(P,h13(2),epsilon3(j,i),VLC_EH3, tau(j));
         end
    end
    S1 = S1 + P11;
    S2 = S2 + P12;
    S3 = S3 + P13;
    
    T1 = T1 + T_trans1;
    T2 = T2 + T_trans2;
    T3 = T3 + T_trans3;
    
end

P1 = S1/loop;
P2 = S2/loop;
P3 = S3/loop;

T1 = T1/loop;
T2 = T2/loop;
T3 = T3/loop;


figure(3);
plot(K,tau(1) - T1(1,:), '--r', K, tau(1) - T2(1,:), ':ob', K, tau(1) - T3(1,:),'-k', K,tau(2) - T1(2,:), '--r',  K, tau(2) - T2(2,:), ':ob', K, tau(2) - T3(2,:),'-k','linewidth',2)
grid on
ylabel('T^{comp}_j (s)')
xlabel('The number of local iterations, K')
legend('User 1', 'User 2', 'User 3')

figure(4);
plot(K,P1(1,:)-P1(1,1), '--r', K, P2(1,:)-P2(1,1), ':ob', K, P3(1,:)-P3(1,1),'-k', K,P1(2,:)-P1(2,1), '--r', K, P2(2,:)-P2(2,1), ':ob', K, P3(2,:)-P3(2,1),'-k','linewidth',2)
grid on
ylabel('The increased amount of transmit IRL power (W)')
xlabel('The number of local iterations, K')
legend('User 1', 'User 2', 'User 3')

