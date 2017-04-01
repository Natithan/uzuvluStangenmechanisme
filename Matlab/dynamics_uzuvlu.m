%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = ...
dynamics_uzuvlu(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
    m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_uzuvlu)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog2_P_x = -X2*cos(phi2)-Y2*cos(phi2+pi/2);
cog2_P_y = -X2*sin(phi2)-Y2*sin(phi2+pi/2);
cog2_Q_x = (r2-X2)*cos(phi2)-Y2*cos(phi2+pi/2);
cog2_Q_y = (r2-X2)*sin(phi2)-Y2*sin(phi2+pi/2);
cog3_Q_x = -Y3*cos(phi3+pi/2)-X3*cos(phi3);
cog3_Q_y = -Y3*sin(phi3+pi/2)-X3*sin(phi3);
cog3_R_x = (r3-X3)*cos(phi3)-Y3*cos(phi3+pi/2);
cog3_R_y = (r3-X3)*sin(phi3)-Y3*sin(phi3+pi/2);
cog4_S_x = -X4*cos(phi4)-Y4*cos(phi4+pi/2);
cog4_S_y = -X4*sin(phi4)-Y4*sin(phi4+pi/2);
cog4_R_x = (r4-X4)*cos(phi4)-Y4*cos(phi4+pi/2);
cog4_R_y = (r4-X4)*sin(phi4)-Y4*sin(phi4+pi/2);

% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];

% 3D model vectors
P_cog2_vec = [-cog2_P_x    -cog2_P_y    zeros(size(phi2))];
Q_cog3_vec = [-cog3_Q_x    -cog3_Q_y    zeros(size(phi2))];
S_cog4_vec = [-cog4_S_x    -cog4_S_y    zeros(size(phi2))];
PQ_vec = [r2*cos(phi2) r2*sin(phi2) zeros(size(phi2))];

% acceleration vectors
acc_2 =       cross(omega2,cross(omega2,P_cog2_vec))+cross(alpha2,P_cog2_vec);
acc_Q =       cross(omega2,cross(omega2,PQ_vec    ))+cross(alpha2,PQ_vec    );
acc_3 = acc_Q+cross(omega3,cross(omega3,Q_cog3_vec))+cross(alpha3,Q_cog3_vec);
acc_4 =       cross(omega4,cross(omega4,S_cog4_vec))+cross(alpha4,S_cog4_vec);
acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);


% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
F_P_x = zeros(size(phi2));
F_P_y = zeros(size(phi2));
F_Q_x = zeros(size(phi2));
F_Q_y = zeros(size(phi2));
F_R_x = zeros(size(phi2));
F_R_y = zeros(size(phi2));
F_S_x = zeros(size(phi2));
F_S_y = zeros(size(phi2));
M_P = zeros(size(phi2));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
  A = [ 1           0            1            0            0            0            0           0           0;
        0           1            0            1            0            0            0           0           0;
        0           0           -1            0           -1            0            0           0           0;
        0           0            0           -1            0           -1            0           0           0;
        0           0            0            0            1            0            1           0           0;
        0           0            0            0            0            1            0           1           0;
       -cog2_P_y(k) cog2_P_x(k) -cog2_Q_y(k)  cog2_Q_x(k)  0            0            0           0           1;
        0           0            cog3_Q_y(k) -cog3_Q_x(k)  cog3_R_y(k) -cog3_R_x(k)  0           0           0;
        0           0            0            0           -cog4_R_y(k)  cog4_R_x(k) -cog4_S_y(k) cog4_S_x(k) 0];
    
  B = [ m2*acc_2x(k);
        m2*acc_2y(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        J2*ddphi2(k);
        J3*ddphi3(k);
        J4*ddphi4(k)];
    
    x = A\B;
    
    % save results
    F_P_x(k) = x(1);
    F_P_y(k) = x(2);
    F_Q_x(k) = x(3);
    F_Q_y(k) = x(4);
    F_R_x(k) = x(5);
    F_R_y(k) = x(6);
    F_S_x(k) = x(7);
    F_S_y(k) = x(8);
    M_P(k)   = x(9);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_P_x,F_P_y),grid
    xlabel('F_P_x [N]')
    ylabel('F_P_y [N]')
    axis tight
    subplot(222)
    plot(F_Q_x,F_Q_y),grid
    xlabel('F_Q_x [N]')
    ylabel('F_Q_y [N]')
    axis tight
    subplot(223)
    plot(F_R_x,F_R_y),grid
    xlabel('F_R_x [N]')
    ylabel('F_R_y [N]')
    axis tight
    subplot(224)
    plot(F_S_x,F_S_y),grid
    xlabel('F_S_x [N]')
    ylabel('F_S_y [N]')
    axis tight
    
    figure
    plot(t,M_P)
    ylabel('M_P [N-m]')
    xlabel('t [s]')
    
end


