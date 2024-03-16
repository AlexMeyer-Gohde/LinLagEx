%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%input_example_mankiw_reis_2007.m
%
%This is the model
%from: 
%Mankiw, N. G., and R. Reis (2007): “Sticky Information in General Equilibrium,” Journal of
%the European Economic Association, 5(2-3), 603–613.
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

%Paramters
v=34.068;
gamma=4.196;
rho_g= 0.93764;%.938;
sigma_g= 0.013929;%.014;
rho_v= 0.62975;%.630;
sigma_v= 1.8194;%1.819;
rho_gamma= .66682;%.667;
sigma_gamma= .18657;%.187;
delta= .18353;%.184;
omega= .19516;%.195;
lambda= .70177;%.702;
beta = 2/3;
psi = 4;
theta = 1;
phi_y = .33;
phi_pi = 1.24;
rho_epsilon = .9179;%.918;
sigma_epsilon = .0116;%.012;
rho_da = .3496;%.350;
sigma_da = .0102;%.010;

%Steady States

%Paramter Abbreviations
e_1w=lambda*beta;
e_1y=lambda*(1-beta);
e_1p=-(beta+(1-lambda)*v*(1-beta));
e_1a=-lambda;
e_1v=-lambda*beta/(v-1);

e_2R=-delta*theta;
e_2yinf=delta;

e_3w=-(psi+gamma*(1-omega));
e_3R=-omega*psi;
e_3l=omega;
e_3p=omega*psi;
e_3yinf=omega*psi/theta;
e_3gam=-omega*psi/(gamma-1);

xi_g=(beta/theta)/(1+1/psi+beta/theta-beta);
xi_v=(beta/(v-1))/(1+1/psi+beta/theta-beta);
xi_gam=(beta/(gamma-1))/(1+1/psi+beta/theta-beta);
xi_a=(1+1/psi)/(1+1/psi+beta/theta-beta);

%Variable names
ENDOGENOUS_VARIABLE_NAMES={'Nominal Wage'
'Output'
'Long Real Interest Rate'
'Labor'
'Nominal Interest Rate'
'Price Level'
'Inflation Rate'
'Natural Output'
'Expected Long Run Natural Output'
'Techonolgy Level'};
EXOGENOUS_VARIABLE_NAMES={'Technology Growth'
'Aggregate Demand'
'Goods Markup'
'Labor Markup'
'Monetary Policy'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0=sum A_iE_t-i(x_t+1)+sum B_iE_t-i(x_t)+sum C_iE_t-i(x_t-1)+sum
% F_iE_t-i(z_t+1)+sumG_iE_t-i(z_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices
%       w   y   R   l   i   p   pi  yn  yninf   a 
A_0=[   0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   1   0   0   0   -1  0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0];
%       w   y   R   l   i   p   pi  yn  yninf   a
B_0=[   e_1w e_1y 0 0   0   e_1p 0  0   0       e_1a
        0   -1  e_2R 0  0   0   0   0   e_2yinf 0
        e_3w 0  e_3R  e_3l 0 e_3p 0 0   e_3yinf 0
        0   -1  0   beta 0  0   0   0   0       1
        0   phi_y 0 0   -1  0   phi_pi -phi_y 0 0
        0   0   0   0   0   1   -1  0   0       0
        0   0   -1  0   1   0   0   0   0       0
        0   0   0   0   0   0   0   -1  0       xi_a
        0   0   0   0   0   0   0   0   0       -1
        0   0   0   0   0   0   0   0   -1      xi_a];
%       w   y   R   l   i   p   pi  yn  yninf   a
C_0=[   0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   -1  0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       0
        0   0   0   0   0   0   0   0   0       1
        0   0   0   0   0   0   0   0   0       0];
%       da  g   v   gam ep
F_0=[	zeros(10,5)];
%       da  g   v   gam ep
G_0=[   0   0   e_1v 0  0
        0   1   0   0   0
        0   0   0   e_3gam 0
        0   0   0   0   0
        0   0   0   0   -1
        0   0   0   0   0
        0   0   0   0   0
        0   xi_g xi_v xi_gam 0
        1   0   0   0   0
        xi_a*rho_da/(1-rho_da) 0 0 0 0];


%name of the counter used in subsequent matrices
it_name='J_1';
%maximum value of counter: either 'infinity' or a scalar (without
%quotation marks)
it_max_value='infinity';


%counter matrices

A_j='[zeros(10,10)]';
%       w   y   R   l   i   p   pi  yn  yninf   a
B_j=['[e_1w*(1-lambda)^J_1 e_1y*(1-lambda)^J_1 0 0 0 lambda*v*(1-beta)*(1-lambda)^J_1 0 0 0 -lambda*(1-lambda)^J_1;'...
     '  0   0   e_2R*(1-delta)^J_1 0 0 0 0 0 e_2yinf*(1-delta)^J_1 0;'...
     'gamma*omega*(1-omega)^J_1 0 e_3R*(1-omega)^J_1 e_3l*(1-omega)^J_1 0 e_3p*(1-omega)^J_1 0 0 e_3yinf*(1-omega)^J_1 0;'...
     'zeros(7,10)]'];
C_j='[zeros(10,10)]';
%       da  g   v   gam ep
G_j=['[  0   0   e_1v*(1-lambda)^J_1 0 0;'...
    ' zeros(1,5);'...
    '   0   0   0   e_3gam*(1-omega)^J_1 0;'...
    'zeros(7,5)]'];
F_j='[zeros(10,5)]';

%limit matrices
%If a counter goes to infinity, please enter the limiting summed matrices
A_inf=A_0;
B_inf=B_0+[e_1w*(1-lambda)/lambda e_1y*(1-lambda)/lambda 0 0 0 v*(1-beta)*(1-lambda) 0 0 0 -(1-lambda);
     0   0   e_2R*(1-delta)/delta 0 0 0 0 0 e_2yinf*(1-delta)/delta 0;
     gamma*(1-omega) 0 e_3R*(1-omega)/omega e_3l*(1-omega)/omega 0 e_3p*(1-omega)/omega 0 0 e_3yinf*(1-omega)/omega 0;
     zeros(7,10)];
C_inf=C_0;
F_inf=F_0;
G_inf=G_0+[  0   0   e_1v*(1-lambda)/lambda 0 0;
    zeros(1,5);
    0   0   0   e_3gam*(1-omega)/omega 0;
    zeros(7,5)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define exogenous process as a VAR(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR-coefficients
%   da  g   v   gam ep
N=[ rho_da 0 0  0   0
    0   rho_g 0 0   0
    0   0   rho_v 0 0
    0   0   0   rho_gamma 0
    0   0   0   0   rho_epsilon];
%Covariance matrix
Omega=[ sigma_da^2 0 0  0   0
        0   sigma_g^2 0 0   0
        0   0   sigma_v^2 0 0
        0   0   0   sigma_gamma^2 0
        0   0   0   0   sigma_epsilon^2];

%Options
PERIOD=4;

linlagex;
