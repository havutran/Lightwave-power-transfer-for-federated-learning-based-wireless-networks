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
K = [4];
tau = [3];

theta = 40*10^(3) ;%Kb -> b

S1=0;S2=0;S3=0;T1=0;T2=0;T3=0;


%----------------------------------------------
%RF channels
Kdb = 8;
dist = 3;
%[g1] = Pathloss_Rician_channels (dist+0.3,1,M, Kdb);
%[g2] = Pathloss_Rician_channels (dist,1,M, Kdb);
%[g3] = Pathloss_Rician_channels (dist-0.3,1,M, Kdb);

g1 = [0.247254151759020 + 0.0515057222263390i;0.229719208822754 + 0.0132816946938692i;0.216181629687725 - 0.0382425794199391i;0.266039936497534 - 0.0360928386299470i];
g2 = [0.297486208883710 + 0.0639872606466984i;0.121674064455455 + 0.0540313306482333i;0.221205232027230 + 7.28649519153710e-05i;0.100542439734184 - 0.00444163453313184i];
g3 = [0.0766562755655168 + 0.00574772609484325i;0.297225272413454 - 0.0682016245633628i;0.0977858285810834 + 0.0295887101162695i;0.0886648640181461 + 0.0486788738422390i];

%----------------------------------------------
%Optimal beamforming
[w1,Gamma_1] = beamforming(g1,g2,g3);
%[w2,Gamma_2] = beamforming(g2,g1,g3);
%[w3,Gamma_3] = beamforming(g3,g2,g1);

alpha = 2*10^(-28);
c = 20;
D = 10*10^(6); %Mb -> bits
fmax = 1.5*10^(9); %GHz -> Hz
fmin = 0.3*10^(9); %GHz -> Hz


%----------------------------------------------
%Solving SubOP1
            a = tau - c*D*K/fmin;
            b = tau - c*D*K/fmax;
            [T_trans1] = goldensectionsearch (a,b,alpha,c,D,theta, Gamma_1, K, tau);
 %----------------------------------------------
%Solving SubOP2           
            %Optimal P1j
            P_max = 120;
            epsilon = fT(T_trans1,alpha,c,D,theta, Gamma_1, K, tau);
            P = bisection(P_max,h11(2),epsilon,VLC_EH1);
            
%----------------------------------------------
%Exhaustive Search           
            test_range = linspace(0.05,2,101);
            
            for i = 1 : 101  
                epsilon1(i) = fT(test_range(i),alpha,c,D,theta, Gamma_1, K, tau);
            end

            for i = 1 : 101  
                P11(i) = bisection(P_max,h11(2),epsilon1(i),VLC_EH1);
            end


figure(5);
plot(test_range, P11,'-k', T_trans1, P,'or', 'linewidth',2)
grid on
ylabel('Transmit IRL Power (W)')
xlabel('Transmission time')
legend('Exhaustive Search', 'Optimal Solution Obtained by Our method')

