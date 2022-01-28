function [CM_cg_matrix, force_coeff, k1, CD0] = LLWing_function(DE_flap,N, BreakCl)

% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% -------------------------------------------------------------------------     

format long;

% -------------------------------------------------------------------------
%% INPUT DATA
% -------------------------------------------------------------------------

rho_inf = 1.225; % air density [kg/m^3]
U_inf = 100/3.6; % plane velocity [m/s]
q_inf = 1/2*rho_inf*U_inf^2;
k = 0.75; % volume correction factor
l_f = 1; % fuselage longitude
r = 0.3; % fuselage radius

inc = 0*pi/180; % fuselage inclination

% Wing planform (assumes planar wing)
%C_T = 0.28; % wing tip chord
%C_R = 1.55; % wing root chord
AR = 21.3 ;   % aspect ratio (AR = span^2/S_alar)
TR = 1/5.55;   % taper ratio (lambda = C_T/C_R)
DE25 = 15; % sweep angle at c/4 (deg) (LAMBDA)
span = 20; % wing span [m]

ETIP = -5.25; % tip twist (deg, negative for washout)

% -------------------------------------------------------------------------
%% AIRFOIL DATA EXTRACTION
% -------------------------------------------------------------------------
[alpha_l0_T, alpha_l0_R, CM0_T, CM0_R, CD0_polin_reg_T, CD0_polin_reg_R] = airfoil_parameters();


% Sections data (uses linear interpolation between root and tip)
A0p = [ alpha_l0_R alpha_l0_T ]; % root and tip section zero-lift angles (deg)
CM0p = [ CM0_R CM0_T ]; % root and tip section free moments
CDP = [ CD0_polin_reg_R(1,3) CD0_polin_reg_R(1,2) CD0_polin_reg_R(1,1) ;   % root section CD0, k1 and k2  (airfoil CD curve)
        CD0_polin_reg_T(1,3) CD0_polin_reg_T(1,2) CD0_polin_reg_T(1,1) ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)
YF_pos = [ 0.0 1.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2;  % flap_chord/chord ratio
FlapCorr = 0.8; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)
%N = 100 ; % number of panels along the span
ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ] ; % angles of attack for analysis (deg)


% -------------------------------------------------------------------------
%% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------

% Wing discretization (lenghts are dimensionless with the wing span)
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr);
% S is the wing surface and mac is the mean aerodynamic chord

% Assembly of the influence coefficients matrix (needed once)
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, but general... and illustrative!)
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

%% Fuselage parameters
[C_MF_adim, C_MF_2_adim, C_D0_adim] = fuselage_parameters(rho_inf, U_inf, k, l_f, r, S, mac, span, ALPHA, inc);

%% Coefficient parametres
[alpha_l0, cl_alpha, mom_slope, x_ac, CM0, CLa, CD0, k1, mean_CM_cg, CM_cg_matrix, CL_trim, CL_stall, v_stall, margin_stab] = coeff_parameters(ALPHA, force_coeff, cl_local, S, span, mac, C_MF_2_adim, C_D0_adim, rho_inf, BreakCl);

close Figure 1 Figure 2 
%Figure 3 Figure 4 Figure 5 Figure 6 Figure 7 Figure 8


end