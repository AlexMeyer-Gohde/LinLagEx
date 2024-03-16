function []=mode_estimates(Pparameter,Ptype,selection_matrix,Data,T,K,model_name,x0,save_name,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%mode_estimates.m
%
%This program:
%               -creates the prior density function using prior means, variances, and types
%               -creates functions transform and inverse_transform to to tramsform the domain of
%                   the parameter to [-Inf, Inf]
%               -finds the approximated posterior mode and inverse Hessian at the mode
%               -calculates and displays some results
%
%Inputs:        Pmean: Prior mean
%               Pvar: Prior variance
%                           Note: for uniform distribution, Pmean and Pvar
%                           contain the bounds of the uniform distribution
%               Ptype: type of prior (norm, beta, gamma, invgamma, or
%                   uniform)
%               Ttype: type of prior used for domain transform (i.e. norm
%                   does not transform the domain, beta has the parameter
%                   restricted to [0,1], etc.
%               selection_matrix: a vector containing an integer at the
%                   i'th entry that denotes which endogenous variable in
%                   the model file corresponds to the i'th observable
%               Data: vectorized data set
%               T: number of observations
%               K: number of observables
%               model_name: a function that executes the model file
%               x0: vector of initial values for numerical optimization
%               save_name: the name of the file to which the results should
%                   be saved.
%
%The code has been adapted from: SIGE_MLE.m,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part of the set of files that accompany the article:               %
%Mankiw, N. Gregory and Ricardo Reis (2007) "Sticky Information in  %
%General Equilibrium," Journal of the European Economic Association,%
%forthcoming.                                                       %
%Last revised: August 30, 2006                                      %
%Written by: Ricardo Reis                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%THIS VERSION: 1.0 March 27, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_name=str2func(model_name);
[log_prior_function]=create_log_priors(Pparameter, Ptype);
[transform inverse_transform]=create_transforms(Pparameter);

if nargin ==9
x0=inverse_transform(x0);

objective_function=@(x) model_loglike(transform(x),selection_matrix,Data,T,K,model_name)-log_prior_function(transform(x));
options = optimset('MaxFunEvals',1000000,'MaxIter', 1000000, 'TolFun', 1e-10, 'TolX', 1e-10, 'HessUpdate','bfgs','LargeScale','off','Display', 'iter');
%[posterior_mode_1 lpd]=fminunc(objective_function,x0,options);

%adding a percent sign to the line above and removing one of those below
%changes the algorithm used for numerical optimization
%[posterior_mode_1,fval,exitflag,output]=fminsearch(objective_function,x0,options);
[lpd,posterior_mode_1,gh,hess,itct,fcount,retcodeh] = csminwel(transform,objective_function,x0,eye(length(x0))*1e-7,[],1e-10,10000);

posterior_mode=transform(posterior_mode_1);
%posterior_mode=transform(x0);

objective_function=@(x) model_loglike(x,selection_matrix,Data,T,K,model_name)-log_prior_function(x);
%options = optimset('MaxFunEvals',1000000,'Display', 'iter', 'TolFun', 1e-10, 'LargeScale', 'off');
%[posterior_mode lpd a b c d hessian]=fmincon(objective_function,x0,[],[],[],[],Pparameter(3,:),Pparameter(4,:),[],options);
%tolerances=[1e-15,1e-2];

lpd=-lpd;
x=posterior_mode';
else
X0=x0;
for j=1:size(X0,1)
disp('iteration')
j
temp=rand(1,length(x0));
x0=X0(j,:);
x0=inverse_transform(x0)+(log(temp)-log(ones(1,length(x0))-temp))./2;
X0(j,:)=transform(x0);
objective_function=@(x) model_loglike(transform(x),selection_matrix,Data,T,K,model_name)-log_prior_function(transform(x));

%options = optimset('MaxFunEvals',500,'MaxIter', 1000, 'TolFun', 1e-10, 'TolX', 1e-10, 'HessUpdate','bfgs','LargeScale','off','Display', 'iter');
%[x0 fval]=fminsearch(objective_function,x0,options);

options = optimset('MaxFunEvals',1000,'MaxIter', 10000, 'TolFun', 1e-5, 'TolX', 1e-5, 'HessUpdate','bfgs','LargeScale','off','Display', 'iter');
[posterior_mode_1 lpd]=fminunc(objective_function,x0,options);

%adding a percent sign to the line above and removing one of those below
%changes the algorithm used for numerical optimization
%[posterior_mode_1,fval,exitflag,output]=fminsearch(objective_function,x0,options);
%[lpd,posterior_mode_1,gh,hess,itct,fcount,retcodeh] = csminwel(transform,objective_function,x0,eye(length(x0))*1e-7,[],1e-4,10000);

posterior_mode(j,:)=[transform(posterior_mode_1) -lpd];
%posterior_mode=transform(x0);
end
posterior_mode
[lpd,row]=max(posterior_mode(:,size(posterior_mode,2)));
x=inverse_transform(posterior_mode(row,1:end-1));
options = optimset('MaxFunEvals',1000000,'MaxIter', 1000000, 'TolFun', 1e-10, 'TolX', 1e-10, 'HessUpdate','bfgs','LargeScale','off','Display', 'iter');
objective_function=@(x) model_loglike(x,selection_matrix,Data,T,K,model_name)-log_prior_function(x);
[posterior_mode_1 lpd]=fminunc(objective_function,x,options);
x=transform(posterior_mode_1)';
lpd=-lpd;
posterior_mode=[posterior_mode;posterior_mode_1 lpd];
end
hess2;
hess=reshape(hessian_mat,[length(x),length(x)]);
%hess=0.5*(hess+hess')
VARCOV = inv(hess)
for j=1:length(x0)
con_interv(j,1:2)=norminv([0.05, 0.95],posterior_mode(j),VARCOV(j,j)^(1/2));
end
laplace_marginal_density=.5*length(x0)*log(2*pi)+.5*log(det(VARCOV))+lpd;
if nargin ==9
save(save_name,'VARCOV','posterior_mode','selection_matrix','Data','T','K','model_name','log_prior_function','Ptype','Pparameter','lpd','laplace_marginal_density')
else
save(save_name,'VARCOV','posterior_mode','selection_matrix','Data','T','K','model_name','log_prior_function','Ptype','Pparameter','X0','lpd','laplace_marginal_density')
posterior_mode=x';
end

display_results=num2str([posterior_mode' diag(VARCOV).^(1/2) con_interv]);
[junk width]=size(display_results);
display_results=[['Mode' blanks(floor((width-12)/3)) 'Std' blanks(floor((width-12)/3)) '5%' blanks(width-12-2*floor((width-12)/3)) '95%']; repmat('-',[1 width]); display_results; repmat('-',[1 width])];
disp(blanks(5)')
disp('Estimation results:')
disp(display_results)
disp(['Log Posterior Density:' blanks(2) num2str(lpd)])
disp(['Laplace Approx of Log Marginal Data Density:' blanks(2) num2str(laplace_marginal_density)])

