%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_uhlig_3.m
%
%This is exampl3.m from Prof. Uhlig's Toolkit
%
% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL3.M calculates through an extensions of Hansens benchmark real business
% cycle model with a simple time-to-decay feature, leading to echo effects.  
% The change to Hansens model consists in assuming that capital
% does not depreciate at all for p periods, and then depreciates 100 percent.  We set p to 4.
% This example was prepared to discuss a paper by Boucekkine, Germain and Licandro
%
% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.
%
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('..\linlagex')

N_bar     = 1.0/3;  % Steady state employment is a third of total time endowment
Z_bar     = 1; % Normalization
rho       = .36; % Capital share
R_bar     = 1.01; % One percent real interest per quarter
eta       = 1.0; % constant of relative risk aversion = 1/(coeff. of intertemporal substitution)
psi       = .95; % autocorrelation of technology shock
sigma_eps = .712; % Standard deviation of technology shock.  Units: Percent.
p_echo    = 4;

%disp('EXAMPLE 3: Hansen benchmark real business cycle model ');
%disp('           with a time-to-decay feature, leading to echo effects.');
%disp(sprintf('          New capital does not depreciate for %2.0f periods',p_echo));
%disp('          but completely depreciates after that');







% Calculating the steady state:

betta   = 1.0/R_bar;  % Discount factor beta
YK_bar  = (1- betta)/((1 - betta^p_echo)*betta*rho); % = Y_bar / K_bar
K_bar   = (YK_bar / Z_bar)^(1.0/(rho-1)) * N_bar;
I_bar   = K_bar / p_echo ;
Y_bar   = YK_bar * K_bar;
C_bar   = Y_bar - I_bar;
Lam_bar = C_bar^(- eta); % Lambda bar, steady state marg. utility of cons.
Mu_bar  = rho*Lam_bar*YK_bar;  % Lagrangian on capital = sum(j=1..p_echo) of lagged investments
A       = Lam_bar * (1 - rho) * Y_bar/N_bar;

%Variable names
ENDOGENOUS_VARIABLE_NAMES={'investment'  % 1
            'investment(t-1)'  % 2
            'investment(t-2)'  % 3
            'investment(t-3)'  % 4
            'E_t[mu(t+2)]'  % 5
            'E_t[mu(t+3)]'  % 6
            'E_t[mu(t+4)]'  % 7
            'consumption'  % 8
            'output'  % 9
            'capital'  % 10
            'labor'  % 11
            'marginal utility'  % 12
            'mu'  % 13
            'E_t[mu(t+1)]'};
EXOGENOUS_VARIABLE_NAMES={'Solow parameter'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices

A_0=[zeros(10,14);[0, 0, 0, 0,    0, 0, 0 0,     0,      0,    0,     0,       -1,        0
       0, 0, 0, 0,    0, 0, 0  0,     0,      0,    0,     0,        0,       -1
       0, 0, 0, 0,   -1, 0, 0    0,     0,      0,    0,     0,        0,        0
       0, 0, 0, 0,    0,-1, 0    0,     0,      0,    0,     0,        0,        0]];
B_0=[-I_bar, 0, 0, 0,             0,        0,        0   -C_bar, Y_bar,      0,    0,     0,        0,        0
     0, 0, 0, 0,             0,        0,        0      0,     0,-p_echo,    0,     0,        0,        0
            0, 0, 0, 0,             0,        0,        0        0,    -1,    rho,(1-rho),   0,        0,        0
            0, 0, 0, 0,             0,        0,        0       0,     1,      0,   -1,     1,        0,        0 
            0, 0, 0, 0,       betta^2,  betta^3,  betta^4       0,     0,      0,0,(-Lam_bar/Mu_bar), 0,    betta 
            0, 0, 0, 0,             0,        0,        0      0,     1,     -1,    0,     1,       -1,        0
            0, 0, 0, 0,             0,        0,        0       eta,     0,      0,    0,     1,        0,        0 
            0, 1, 0, 0,             0,        0,        0        0,     0,      0,    0,     0,        0,        0
            0, 0, 1, 0,             0,        0,        0       0,     0,      0,    0,     0,        0,        0
            0, 0, 0, 1,             0,        0,        0   0,     0,      0,    0,     0,        0,        0
        0, 0, 0, 0,    0, 0, 0      0,     0,      0,    0,     0,        0,        1
       0, 0, 0, 0,    1, 0, 0   0,     0,      0,    0,     0,        0,        0
       0, 0, 0, 0,    0, 1, 0  0,     0,      0,    0,     0,        0,        0
       0, 0, 0, 0,    0, 0, 1    0,     0,      0,    0,     0,        0,        0];

C_0=[[0, 0, 0, 0,             0,        0,        0   
            1, 1, 1, 1,             0,        0,        0   
            0, 0, 0, 0,             0,        0,        0   
            0, 0, 0, 0,             0,        0,        0   
            0, 0, 0, 0,             0,        0,        0   
            0, 0, 0, 0,             0,        0,        0   
            0, 0, 0, 0,             0,        0,        0  
           -1, 0, 0, 0,             0,        0,        0  
            0,-1, 0. 0,             0,        0,        0   
            0, 0,-1, 0,             0,        0,        0 
        0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0
       0, 0, 0, 0,    0, 0, 0],zeros(14,7)];
F_0=[zeros(14,1)];
G_0=[ 0
       0
       1
       0
       0
       0
       0
       0
       0
       0
zeros(4,1)];


%name of the counter used in subsequent matrices
%it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value=0;


%counter matrices



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
N=[psi];
%Covariance matrix
Omega=[sigma_eps^2];

%Options
PERIOD=4;


linlagex;
