function [loglik, varargout]=model_loglike(parameters,Q,Data,T,K,model_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%model_loglike.m
%
%This program:
%               prepares the block Levinson algorithm for computing the
%               loglikelihood by calling on the model to get the recusive
%               solution and then, using spectral methods, generate the
%               sequence of autocovariance matrices
%
%Input:         parameters: a vector of parameters to be passed to the
%                   model
%               Q: a cell array such that
%                   Observables=(Q{1}-Q{2}L)*Endogenous is stationary,
%                   where L is the lag operator and Q{1} and Q{2} are
%                   matrices of appropriate dimensions (K x num_endog)
%                   - if the transform is a simple linear transform, let
%                   Q{1} be such that Observables=Q{1}*Endogenous and
%                   define Q{2}=[];
%               Data: vectorized data set
%               T: number of observations
%               K: number of observables
%               model_name: a function that executes the model file
%
%Output:        loglik: the value of the log-likelihood function with the
%                   given parameters and data set 
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


existence_uniqueness=1;
trisolve=1;
model_name();
if isnumeric(existence_uniqueness)==1
if isempty(Q{2})==0
Q_1=Q{1}(:,:);
Q_2=Q{2}(:,:);

R_1=Q_1;    
R_2=Q_2-R_1*ALPHA_ZS;

grid_size=2*T;
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;


s_Y=cell(1,grid_size);
s_yy=zeros(K*K,grid_size);

THETA=reshape(MA_VECTOR(1:num_endog*max_it,:), [ num_endog max_it num_exog ] );
THETA=permute(THETA, [1 3 2]);
THETA=reshape(THETA,[num_endog num_exog*(max_it)]);
vec_THETA=reshape(THETA,[num_endog*num_exog,numel(THETA)/(num_endog*num_exog)]);
if max_it>0
THETA_I_MINUS_1=MA_VECTOR((max_it-1)*num_endog+1:(max_it)*num_endog,:);
else
THETA_I_MINUS_1=zeros(num_endog,num_exog);
end

dd_1=repmat(grid_points', [1 max_it]);%
dd_2=repmat([0:max_it-1], [length(grid_points) 1]);%
dd=dd_1.*dd_2;
E_plus=exp( ((-1)^(1/2))*grid_points);
EE_PLUS=exp( ((-1)^(1/2))*dd);
%%%%chol_Omega=chol(Omega);
for n = 1 : grid_size
e_plus=E_plus(n);
e_minus=conj(e_plus);
THETA_E_MINUS=vec_THETA*EE_PLUS(n,:)';
THETA_E_MINUS=(Q_1-Q_2*e_minus)*reshape(THETA_E_MINUS,[num_endog num_exog]);
Temp=(N^(max_it))/(eye(num_exog)-N*e_minus);
Temp=(R_1-R_2*e_minus)*BETA_ZS*Temp;
Temp=Temp+R_1*ALPHA_ZS*THETA_I_MINUS_1;
%Temp=(eye(num_endog)-ALPHA_ZS*e_minus)\Temp;
Temp=Temp*(e_minus^(max_it));
Temp=THETA_E_MINUS+Temp;
%%%Temp=Temp*chol_Omega;
%Temp=selection_matrix*Temp;
Temp=((2*pi)^(-1))*Temp*Omega*Temp';
s_Y{n}=Temp(:);
end
else
Q_1=Q{1}(:,:);

grid_size=2*T;
grid_points=0:2*pi/grid_size:2*pi-2*pi/grid_size;


s_Y=cell(1,grid_size);
s_yy=zeros(K*K,grid_size);

THETA=reshape(MA_VECTOR(1:num_endog*max_it,:), [ num_endog max_it num_exog ] );
THETA=permute(THETA, [1 3 2]);
THETA=reshape(THETA,[num_endog num_exog*(max_it)]);
vec_THETA=reshape(THETA,[num_endog*num_exog,numel(THETA)/(num_endog*num_exog)]);
if max_it>0
THETA_I_MINUS_1=MA_VECTOR((max_it-1)*num_endog+1:(max_it)*num_endog,:);
else
THETA_I_MINUS_1=zeros(num_endog,num_exog);
end

dd_1=repmat(grid_points', [1 max_it]);%
dd_2=repmat([0:max_it-1], [length(grid_points) 1]);%
dd=dd_1.*dd_2;
E_plus=exp( ((-1)^(1/2))*grid_points);
EE_PLUS=exp( ((-1)^(1/2))*dd);
for n = 1 : grid_size
e_plus=E_plus(n);
e_minus=conj(e_plus);
THETA_E_MINUS=vec_THETA*EE_PLUS(n,:)';
THETA_E_MINUS=reshape(THETA_E_MINUS,[num_endog num_exog]);
Temp=(N^(max_it))/(eye(num_exog)-N*e_minus);
Temp=BETA_ZS*Temp;
Temp=Temp+ALPHA_ZS*THETA_I_MINUS_1;
Temp=(eye(num_endog)-ALPHA_ZS*e_minus)\Temp;
Temp=Temp*(e_minus^(max_it));
Temp=THETA_E_MINUS+Temp;
Temp=Q_1*Temp;
Temp=((2*pi)^(-1))*Temp*Omega*Temp';
s_Y{n}=Temp(:);
end




end

s_yy=cell2mat(s_Y);

s_yy=ifft(conj(s_yy'),'symmetric')*2*pi;
s_yy=reshape(s_yy', [K K*grid_size]);



loglik = block_levinson_loglik(Data, s_yy, 1); 
loglik=-loglik;
else
loglik = 10^50;
end
varargout{1}=max_it;