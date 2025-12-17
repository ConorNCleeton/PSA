%% Script file to run a single cycle of the PSA Simulation code
clc; clear; close all;

format long
load('ParamsNew')
addpath('CycleSteps')
N = 10 ;                        % Number of finite volume elements
type = 'ProcessEvaluation' ;    % Purity / Recovery, but 
                                % not necessary in this code (not
                                % optimising cycle parameters).

%% ************** INDIVIDUAL CYCLE ********************************
% Load in the required material data for zeolite 13X
IsothermParams = [3.36, 2.6, 5.1e-10, 9.76e-12,-34600,-34500,3.36,2.6,2.59e-10,2.59e-10,-19120,-19120,0];
q_s_b      = [IsothermParams(1),  IsothermParams(7)]   ;   % Saturation loading on site b [mol/kg]
q_s_d      = [IsothermParams(2),  IsothermParams(8)]   ;   % Saturation loading on site d [mol/kg]
b          = [IsothermParams(3),  IsothermParams(9)]   ;   % Pre-exponential factor for site b [Pa-1]
d          = [IsothermParams(4),  IsothermParams(10)]  ;   % Pre-exponential factor for site d [Pa-1]
deltaU_b   = [IsothermParams(5),  IsothermParams(11)]  ;   % Heat of adsorption for site b [J/mol]
deltaU_d   = [IsothermParams(6),  IsothermParams(12)]  ;   % Heat of adsorption for site d [J/mol]
IsothermParameters = [q_s_b, q_s_d, b, d, deltaU_b, deltaU_d, IsothermParams(13)] ;

density            = 1730;
cp_s               = 1070;
material_property  = [density, cp_s];

% fixed cycle configuration (other process params in 
% ProcessInputParameters.m file
P_H = 1.0e5;
t_feed = 1000;
alpha_LR = 0.069;
v_feed = 1.27;
epsilon = 0.37;
r_p = 1e-3;
e_p = 0.35;

x = [P_H, v_feed, t_feed, alpha_LR, epsilon, e_p, r_p];

% Convert x into the process variables which are fed into the
% simulation code. The values of x (which were the optimisable cycle
% parameters in the optimisation code) are obtained from a Process
% Optimisation of the Purity Recovery.
% P_0         Adsorption pressure [Pa]
% ndot_0      Inlet molar flux [mol/s/m^2]
% t_ads       Time of adsorption step [s]
% alpha       Light product reflux ratio [-]
% epsilon     Void fraction of the bed [-]
% r_p         Radius of the pellets [m]
% e_p         pellet porosity [-]
% ro_crys     crystal density [kg / m3]
% c_ps        crystal specific heat capacity [J / kg / K]

process_variables = [x(1), x(1)*x(2)/8.314/313.15,  x(3),  x(4),  x(5), x(6),x(7)] ;
y_0              = 0.15             ; 
q                = Isotherm(y_0, P_H, 313.15, IsothermParameters) ;
x0               = zeros(5*N+10,1) ;
x0(1:N+2)        = P_H/P_H         ;
x0(N+3)          = y_0             ;
x0(N+4:2*N+4)    = y_0             ;
x0(2*N+5:3*N+6)  = q(1)/(5.84*density*(1-e_p))      ;
x0(3*N+7:4*N+8)  = q(2)/(5.84*density*(1-e_p))      ;
x0(4*N+9)        = 1               ;
x0(4*N+10:5*N+10)= 313.15/313.15   ;

[objs,~,a,b,c,d,e] = PSACycle_1iter(process_variables, IsothermParams, material_property, type, N,x0) ;
purity = -objs(1);
recovery = -objs(2);
 
%% Post-processing
%      a-e: These are the non-dimensional state variables for the five steps: 
%      a:CoCPressurization, b:Adsorption, c:Heavy Reflux, d:CnCDepressurization 
%      and e:Light Reflux.

%% Dimensional pressures 
% (multiply by P_ads [x(1)])
figure('Name','State Variables','Color','w')
% subplot(2,2,1)
% plot((0:N+1)/(N+1), a(end,1:N+2).*x(1),'-or')
% hold on
% plot((0:N+1)/(N+1), b(end,1:N+2).*x(1),'-ob')
% plot((0:N+1)/(N+1), c(end,1:N+2).*x(1),'-ok')
% plot((0:N+1)/(N+1), d(end,1:N+2).*x(1),'-om')
% plot((0:N+1)/(N+1), e(end,1:N+2).*x(1),'-oc')
% legend('Pres','Ad','HR','CnCDepres','LR','Location','best')
% xlabel('Z / -')
% ylabel('Pressure / Pa')

% Dimensional temperature 
% 4*N+9:5*N+10
% multiply the column temperature by the feed temperature, T_0
T_0 = 313.15;
subplot(2,2,1)
plot((0:N+1)/(N+1), a(end,4*N+9:5*N+10).*T_0,'-or')
hold on
plot((0:N+1)/(N+1), b(end,4*N+9:5*N+10).*T_0,'-ob')
plot((0:N+1)/(N+1), c(end,4*N+9:5*N+10).*T_0,'-ok')
plot((0:N+1)/(N+1), d(end,4*N+9:5*N+10).*T_0,'-om')
plot((0:N+1)/(N+1), e(end,4*N+9:5*N+10).*T_0,'-oc')
legend('Pres','Ad','HR','CnCDepres','LR','Location','best')
xlabel('Z / -')
ylabel('Pressure / Pa')

% CO2 gas mole fraction
% N+3:2*N+4
subplot(2,2,2)
plot((0:N+1)/(N+1), a(end,(N+3):(2*N+4)),'-or')  
hold on
plot((0:N+1)/(N+1), b(end,(N+3):(2*N+4)),'-ob')
plot((0:N+1)/(N+1), c(end,(N+3):(2*N+4)),'-ok')
plot((0:N+1)/(N+1),d(end,(N+3):(2*N+4)),'-om')
plot((0:N+1)/(N+1),e(end,(N+3):(2*N+4)),'-oc')
legend('Pres','Ad','HR','CnCDepres','LR','Location','best')
xlabel('Z / -')
ylabel('y_{CO_2}')


% Dimensional CO2 molar loadings 
% 2*N+5:3*N+6
% multiply the dimensionless molar loading by the molar loading scaling
% factor, q_s

ro_s = density;                % density of adsorbent [kg/m3]
q_s  = 5.84                ;   % Molar loading scaling factor [mol/kg]
q_s0 = ro_s*q_s;               % Molar loading scaling factor [mol/m3]

subplot(2,2,3)
plot((0:N+1)/(N+1),a(end,(2*N+5):(3*N+6)).*q_s0,'-or')
hold on 
plot((0:N+1)/(N+1),b(end,(2*N+5):(3*N+6)).*q_s0,'-ob')
plot((0:N+1)/(N+1),c(end,(2*N+5):(3*N+6)).*q_s0,'-ok')
plot((0:N+1)/(N+1),d(end,(2*N+5):(3*N+6)).*q_s0,'-om')
plot((0:N+1)/(N+1),e(end,(2*N+5):(3*N+6)).*q_s0,'-oc')
legend('Pres','Ad','HR','CnCDepres','LR','Location','best')
xlabel('Z / -')
ylabel('Molar loading of CO_2 in bed / mol m^{-3}')


% Dimensional N2 molar loadings
% 3*N+7:4*N+8
subplot(2,2,4)
plot((0:N+1)/(N+1),a(end,(3*N+7):(4*N+8)).*q_s0,'-or')
hold on 
plot((0:N+1)/(N+1),b(end,(3*N+7):(4*N+8)).*q_s0,'-ob')
plot((0:N+1)/(N+1),c(end,(3*N+7):(4*N+8)).*q_s0,'-ok')
plot((0:N+1)/(N+1),d(end,(3*N+7):(4*N+8)).*q_s0,'-om')
plot((0:N+1)/(N+1),e(end,(3*N+7):(4*N+8)).*q_s0,'-oc')
legend('Pres','Ad','HR','CnCDepres','LR','Location','best')
xlabel('Z / -')
ylabel('Molar loading of N_2 in bed / mol m^{-3}')


