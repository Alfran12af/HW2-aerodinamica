function [alpha_l0_T, alpha_l0_R, CM0_T, CM0_R, CD0_polin_reg_T, CD0_polin_reg_R] = airfoil_parameters()
% Re=600000 is used for wing tip
% Re=1000000 is used for wing root


load('airfoil_data_workspace.mat') % open imported airfoil (CH60) data from csv file


%% Extract csv data
A = table2array(mh60pol1);
Re = [100000; 200000; 400000; 600000; 1000000];


%% Cl as a function of alpha
for i = 1:5 % compute every reynolds number
    cl_matrix(i,:) = A(i,2:8:end); % create cl matrix
    alpha_matrix(i,:) = A(i,1:8:end); %create alpha matrix
end

% linear regression for Re=600000
lin_reg_T = polyfit(A(4,1:8:176), A(4,2:8:176), 1); % wing tip regression line
Cl_alpha_T = lin_reg_T(1,1);
alpha_l0_T = -lin_reg_T(1,2)/lin_reg_T(1,1);

% and Re=1000000
lin_reg_R = polyfit(A(5,1:8:176), A(5,2:8:176), 1); % wing root regression line
Cl_alpha_R = lin_reg_R(1,1);
alpha_l0_R = -lin_reg_R(1,2)/lin_reg_R(1,1);


%% Cd as a function of Cl
for i = 1:5 % compute every reynolds number
    cd_matrix(i,:) = A(i,3:8:end); % create cd matrix
end

% curva regresión para el cd
% Re=600000
CD0_polin_reg_T = polyfit(A(4,2:8:176), A(4,3:8:176), 2); % wing tip
CD0_T = CD0_polin_reg_T(1,3);

% Re=1000000
CD0_polin_reg_R = polyfit(A(5,2:8:176), A(5,3:8:176), 2); % wing root
CD0_R = CD0_polin_reg_R(1,3);


%% Moment coefficient
% the given moment coefficient is refered to the aerodynamic center

for i = 1:5 % compute every reynolds number
    cm_ac_matrix(i,:) = A(i,4:8:end); % create cm matrix
end

% compute the moment coefficient from ac for every Reynolds number
cm_ac_matrix_T = cm_ac_matrix(4,:); % tip free moment coefficient
cm_ac_matrix_R = cm_ac_matrix(5,:); % root free moment coefficient

% compute the mean moment coefficient from ac for the wing tip
CM0_T = mean(cm_ac_matrix_T(1,1:22));

% compute the mean moment coefficient from ac for the wing root
CM0_R = mean(cm_ac_matrix_R(1,1:22));

end