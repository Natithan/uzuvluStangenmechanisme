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
            kinematics_uzuvlu(r1_24,r1_46,h,r2,r3,r4_13,r4_15,r4_35,r5,r6,r7,phi1_14,phi1_46,phi2,dphi2,ddphi2,phi3_init,phi4_13_init,phi5_init,phi6_init,phi7_init,X_init,t,fig_kin_uzuvlu)

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
ddphi4_15 = zeros(size(t));
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
    [X, fval, exitflag]=fsolve('loop_closure_eqs',[phi3_init phi4_13_init,phi5_init,phi6_init,phi7_init,X_init ]',optim_options,phi2(k),r1_24,r1_46,r2,r3,r4_13,r4_15,r4_35,r5,r6,r7,h,phi1_14,phi1_46);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi3(k)=X(1);
    phi4_13(k)=X(2);
    phi5(k)=X(3);
    phi6(k)=X(4);
    phi7(k)=X(5);
    X(k)=X(6);
    
    phi4_15 = zeros(size(t));
    phi4_15(k)=phi4_13(k) - deltaphi13_15;
    
    % *** velocity analysis ***
    A = [-r3*sin(phi3(k)),  r4_13*sin(phi4_13(k)), 0, 0, 0, 0;
         r3*cos(phi3(k)), -r4_13*cos(phi4_13(k)), 0, 0, 0, 0;
         0,0, -r5*sin(phi5(k)),  r6*sin(phi6(k)), 0, 0;
        0,0, r5*cos(phi5(k)), -r6*cos(phi6(k)),0, 0;
        0, 0, 0, 0, r7*sin(phi7(k)),  1;
         0, 0, 0, 0, -r7*cos(phi7(k)), 0];
    B = [ r2*sin(phi2(k))*dphi2(k);
         -r2*cos(phi2(k))*dphi2(k);
        r4_15*sin(phi4_15(k))*dphi4_13(k);%equals dphi4_15(k)
         -r4_15*cos(phi4_15(k))*dphi4_13(k);%equals dphi4_15(k)
     r6*sin(phi6(k))*dphi6(k);
         -r6*cos(phi6(k))*dphi6(k)];
    
    x = A\B;
    
    % save results
    dphi3(k) = x(1);
    dphi4_13(k) = x(2); %equals dphi4_15(k)
    
    dphi5(k) = x(3);
    dphi6(k) = x(4);
    
    dphi7(k) = x(5);
    dX(k) = x(4);
    
    
    % *** acceleration analysis ***
    
    A = [-r3*sin(phi3(k)),  r4_13*sin(phi4_13(k)), 0, 0, 0, 0;
         r3*cos(phi3(k)), -r4_13*cos(phi4_13(k)), 0, 0, 0, 0;
         0,0, -r5*sin(phi5(k)),  r6*sin(phi6(k)), 0, 0;
        0,0, r5*cos(phi5(k)), -r6*cos(phi6(k)),0, 0;
        0, 0, 0, 0, r7*sin(phi7(k)),  1;
         0, 0, 0, 0, -r7*cos(phi7(k)), 0];
    B = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4_13*cos(phi4_13(k))*dphi4_13(k)^2;
         r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4_13*sin(phi4_13(k))*dphi4_13(k)^2;
         
         r4_15*cos(phi4_15(k))*dphi4_13(k)^2+r4_15*sin(phi4_15(k))*ddphi4_15(k)+r5*cos(phi5(k))*dphi5(k)^2-r6*cos(phi6(k))*dphi6(k)^2;
         r4_15*sin(phi4_15(k))*dphi4_13(k)^2-r4_15*cos(phi4_15(k))*ddphi4_15(k)+r5*sin(phi5(k))*dphi5(k)^2-r6*sin(phi6(k))*dphi6(k)^2;
         
         r6*cos(phi6(k))*dphi6(k)^2+r6*sin(phi6(k))*ddphi6(k)-r7*cos(phi7(k))*dphi7(k)^2-0;
         r6*sin(phi6(k))*dphi6(k)^2-r6*cos(phi6(k))*ddphi6(k)-r7*sin(phi7(k))*dphi7(k)^2-0];
    
    x = A\B;
    % save results
    ddphi3(k) = x(1);
    ddphi4_15(k) = x(2); %equals dphi4_15(k)
    
    ddphi5(k) = x(3);
    ddphi6(k) = x(4);
    
    ddphi7(k) = x(5);
    ddX(k) = x(6);
    
    
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
H = r1_24*exp(j*phi1_14);
% point F = fixed
F = r1_24*exp(j*phi1_14) + r1_46*exp(j*phi1_46);
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
    G = F + [X(index),h];
    
    loop1 = [A B C H];
    loop2 = [H D E F];
    loop3 = [F E G];
    
    figure(10)
    clf %Clear current figure window
    hold on
    plot(real(loop1),imag(loop1),'-o',real(loop2),imag(loop2),'-o',real(loop3),imag(loop3),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_uzuvlu
    
    %plot assembly at a certain timestep 
    index = 10; %select 1st timestep
    B = A + r2 * exp(1i*phi2(index));
    C = B + r3 * exp(1i*phi3(index));
    deltaphi13_15 = acos(1/(2*r4_13*r4_15)*(r4_13^2+r4_15^2-r4_35^2));
    phi4_15=phi4_13 - deltaphi13_15;
    D = H + r4_15 * exp(1i*phi4_15(index));
    E = F + r6 * exp(1i*phi6(index));
    G = E + r7 * exp(1i*phi7(index));
    
    figure
    assembly1=[A, B, C, D, E, G];
    assembly2=[C, D, H, C];
    assembly3=[E, F];
    plot(real(assembly1),imag(assembly1),'ro-',real(assembly2),imag(assembly2),'ro-',real(assembly3),imag(assembly3),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis(movie_axis)
    
    figure
    subplot(311)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(312)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(313)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
end



