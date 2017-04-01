%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Analyse van een uzuvlu mechanisme
%
% Nathan Cornille <nathan.cornille@student.kuleuven.be>
% Sander Mergan <sander.mergan@student.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs(phi_init,phi2,r1_24,r1_46,r2,r3,r4_13,r4_15,r4_35,r5,r6,r7,h,phi1_14,phi1_46)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles
phi3=phi_init(1);
phi4_13=phi_init(2);
phi5=phi_init(3);
phi6=phi_init(4);
phi7=phi_init(5);
X=phi_init(6);

deltaphi13_15 = acos(1/(2*r4_13*r4_15)*(r4_13^2+r4_15^2-r4_35^2));
phi4_15=phi4_13 - deltaphi13_15;


% loop closure equations:
F(1)=r2*cos(phi2)+r3*cos(phi3)-r4_13*cos(phi4_13)-r1_24*cos(phi1_14); %L1X
F(2)=r2*sin(phi2)+r3*sin(phi3)-r4_13*sin(phi4_13)-r1_24*sin(phi1_14); %L1Y

F(3)=r4_15*cos(phi4_15)+r5*cos(phi5)-r6*cos(phi6)-r1_46*cos(phi1_46); %L2X
F(4)=r4_15*sin(phi4_15)+r5*sin(phi5)-r6*sin(phi6)-r1_46*sin(phi1_46); %L2Y

F(5)=r6*cos(phi6)-r7*cos(phi7)-X; %L3X
F(6)=r6*sin(phi6)-r7*sin(phi7)-h; %L3Y

