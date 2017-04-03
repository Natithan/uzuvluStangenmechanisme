%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Analyse van een uzuvlumechanisme
%
%
% Nathan Cornille <nathan.cornille@student.kuleuven.be>
% Sander Mergan <sander.mergan@student.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi3,phi4_13,phi5,phi6,phi7,X,dphi3,dphi4_13,dphi5,dphi6,dphi7,dX,ddphi3,ddphi4_13,ddphi5,ddphi6,ddphi7,ddX] = ... 
            kinematics_uzuvlu(r1_24,r1_46,h,r2,r3,r4_13,r4_15,r4_35,r5,r6,r7,phi1_24,phi1_46,phi2,dphi2,ddphi2,phi3_init,phi4_13_init,phi5_init,phi6_init,phi7_init,X_init,t,fig_kin_uzuvlu)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.


phi3 = zeros(size(t));
phi4_13 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
X = zeros(size(t));
dphi3 = zeros(size(t));
dphi4_13 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dX = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi4_13 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddX = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps

deltaphi13_15 = acos(1/(2*r4_13*r4_15)*(r4_13^2+r4_15^2-r4_35^2));

for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument r1_24 ... phi1_68: constants
    % return value x: solution for the unknown angles phi3 and others
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi3_init phi4_13_init,phi5_init,phi6_init,phi7_init,X_init ]',optim_options,phi2(k),r1_24,r1_46,r2,r3,r4_13,r4_15,r4_35,r5,r6,r7,h,phi1_24,phi1_46);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi3(k)=x(1);
    phi4_13(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi7(k)=x(5);
    X(k)=x(6);
    
    phi4_15 = zeros(size(t));
    phi4_15(k)=phi4_13(k) - deltaphi13_15;
    
    % *** velocity analysis ***
    A = [-r3*sin(phi3(k)),  r4_13*sin(phi4_13(k)),  0,                0,                0,                0;
          r3*cos(phi3(k)), -r4_13*cos(phi4_13(k)),  0,                0,                0,                0;
          0,               -r4_15*sin(phi4_15(k)), -r5*sin(phi5(k)),  r6*sin(phi6(k)),  0,                0;
          0,                r4_15*cos(phi4_15(k)),  r5*cos(phi5(k)), -r6*cos(phi6(k)),  0,                0;
          0,                0,                      0,                r6*sin(phi6(k)), -r7*sin(phi7(k)),  1;
          0,                0,                      0,                r6*cos(phi6(k)), -r7*cos(phi7(k)),  0];
    B = [ r2*sin(phi2(k))*dphi2(k);
         -r2*cos(phi2(k))*dphi2(k);
          0;
          0;
          0;
          0];
    
    q = A\B;
    
    % save results
    dphi3(k) = q(1);
    dphi4_13(k) = q(2); %equals dphi4_15(k)
    
    dphi5(k) = q(3);
    dphi6(k) = q(4);
    
    dphi7(k) = q(5);
    dX(k) = q(6);
    
    
    % *** acceleration analysis ***
    
    A = [-r3*sin(phi3(k)),  r4_13*sin(phi4_13(k)),  0,                0,                0,                0;
          r3*cos(phi3(k)), -r4_13*cos(phi4_13(k)),  0,                0,                0,                0;
          0,               -r4_15*sin(phi4_15(k)), -r5*sin(phi5(k)),  r6*sin(phi6(k)),  0,                0;
          0,                r4_15*cos(phi4_15(k)),  r5*cos(phi5(k)), -r6*cos(phi6(k)),  0,                0;
          0,                0,                      0,                r6*sin(phi6(k)), -r7*sin(phi7(k)),  1;
          0,                0,                      0,                r6*cos(phi6(k)), -r7*cos(phi7(k)),  0];
    B = [ r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4_13*cos(phi4_13(k))*dphi4_13(k)^2;
          r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4_13*sin(phi4_13(k))*dphi4_13(k)^2;
         
          r4_15*cos(phi4_15(k))*dphi4_13(k)^2+r5*cos(phi5(k))*dphi5(k)^2-r6*cos(phi6(k))*dphi6(k)^2;
          r4_15*sin(phi4_15(k))*dphi4_13(k)^2+r5*sin(phi5(k))*dphi5(k)^2-r6*sin(phi6(k))*dphi6(k)^2;
         
         -r6*cos(phi6(k))*dphi6(k)^2+r7*cos(phi7(k))*dphi7(k)^2-0;
          r6*sin(phi6(k))*dphi6(k)^2-r7*sin(phi7(k))*dphi7(k)^2-0];
    
    q = A\B;
    % save results
    ddphi3(k) = q(1);
    ddphi4_13(k) = q(2); %equals dphi4_15(k)
    
    ddphi5(k) = q(3);
    ddphi6(k) = q(4);
    
    ddphi7(k) = q(5);
    ddX(k) = q(6);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_13_init = phi4_13(k)+Ts*dphi4_13(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    X_init = X(k)+Ts*dX(k);
    
    
end % loop over positions



% *** create movie ***

% point A = fixed
A = 0;
% point H = fixed
H = r1_24*exp(1j*phi1_24);
% point F = fixed
F = r1_24*exp(1j*phi1_24) + r1_46*exp(1j*phi1_46);
% define which positions we want as frames in our movie
frames = 100;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -r2*2;
y_bottom = -r2*2;
x_right = 10*r2*1.5;
y_top = 5*r2*1.5;

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    B = A + r2 * exp(1i*phi2(index));
    C = B + r3 * exp(1i*phi3(index));
    deltaphi13_15 = acos(1/(2*r4_13*r4_15)*(r4_13^2+r4_15^2-r4_35^2));
    phi4_15=phi4_13 - deltaphi13_15;
    D = H + r4_15 * exp(1i*phi4_15(index));
    E = F + r6 * exp(1i*phi6(index));
    G = F + X(index) + h*1i;
    
    loop1 = [A B C D E G];
    loop2 = [C D H C];
    loop3 = [F E];
    
    figure(10)
    clf %Clear current figure window
    hold on
    plot(real(loop1),imag(loop1),'-o',real(loop2),imag(loop2),'-o',real(loop3),imag(loop3),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save uzuvlu_movie Movie
close(10)


% *** plot figures ***

if fig_kin_uzuvlu
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    B = A + r2 * exp(1i*phi2(index));
    C = B + r3 * exp(1i*phi3(index));
    deltaphi13_15 = acos(1/(2*r4_13*r4_15)*(r4_13^2+r4_15^2-r4_35^2));
    phi4_15=phi4_13 - deltaphi13_15;
    D = H + r4_15 * exp(1i*phi4_15(index));
    E = F + r6 * exp(1i*phi6(index));
    G = F + X(index)+1i*h;
    
    figure
    assembly1=[A, B, C, D, E, G];
    assembly2=[C, D, H, C];
    assembly3=[F, E];
    plot(real(assembly1),imag(assembly1),'ro-',real(assembly2),imag(assembly2),'ro-',real(assembly3),imag(assembly3),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis(movie_axes)
    
%     temp = abs(r2*exp(1i*phi2)+r3*exp(1i*phi3) - (r1_24*exp(1i*phi1_24)+r4_13*exp(1i*phi4_13)))
%     figure
%     plot(t,temp)
%     ylabel('verschil op C via ABC en via AHC [rad]')
%     title('Controle punt C')
%     xlabel('t [s]')
    
    figure
    subplot(421)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(422)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(423)
    plot(t,phi4_13)
    ylabel('\phi_4_13 [rad]')
    subplot(424)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    subplot(425)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    subplot(426)
    plot(t,phi7)
    ylabel('\phi_7 [rad]')
    subplot(428)
    plot(t,X)
    ylabel('X [rad]')
    xlabel('t [s]')
    
    figure
    subplot(421)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(422)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(423)
    plot(t,dphi4_13)
    ylabel('d\phi_4_13 [rad/s]')
    subplot(424)
    plot(t,dphi5)
    ylabel('d\phi_5 [rad/s]')
    subplot(425)
    plot(t,dphi6)
    ylabel('d\phi_6 [rad/s]')
    subplot(426)
    plot(t,dphi7)
    ylabel('d\phi_7 [rad/s]')
    subplot(428)
    plot(t,dX)
    ylabel('dX [rad]')
    xlabel('t [s]')
    
    figure
    subplot(421)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(422)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(423)
    plot(t,ddphi4_13)
    ylabel('dd\phi_4_13 [rad/s^2]')
    subplot(424)
    plot(t,ddphi5)
    ylabel('dd\phi_5 [rad/s^2]')
    subplot(425)
    plot(t,ddphi6)
    ylabel('dd\phi_6 [rad/s^2]')
    subplot(426)
    plot(t,ddphi7)
    ylabel('dd\phi_7 [rad/s^2]')
    subplot(428)
    plot(t,ddX)
    ylabel('ddX [rad]')
    xlabel('t [s]')
    
end



