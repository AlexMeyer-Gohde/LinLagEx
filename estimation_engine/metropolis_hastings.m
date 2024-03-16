function [results,draws,results_densities,lagged_expectations]=metropolis_hastings(mode_estimate_data, draw_type, number_of_draws, burn_in, overdispersion,save_name,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%metropolis_hastings.m
%
%This program draws from the unknown posterior distribution using
%MCMC methods.
%
%Input:     mode_estimate_data: file where the results of the numerical maximization of the
%               posterior mode are located
%           draw_type: 0 if the chain is to be started at the posterior
%               mode, 1 if the chain should be started from a random point
%               draw from the distribution MVNormal (poseterior_mode,
%               inv(hessian)_at_posterior_mode)
%           number_of_draws: the number of draws
%           burn_in: the burn in or number of initial draws to be discarded
%           overdispersion: a parameter to disperse the covariance matrix
%               for draws, change this to get sensible acceptance rates
%           save_name: the name of the file to which the results should be
%               saved
%
%
%The code has been adapted from: SIGE_Bayes.m, SIGE_metro.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part of the set of files that accompany the article:               %
%Mankiw, N. Gregory and Ricardo Reis (2007) "Sticky Information in  %
%General Equilibrium," Journal of the European Economic Association,%
%forthcoming.                                                       %
%Last revised: August 30, 2006                                      %
%Written by: Ricardo Reis                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

load(mode_estimate_data);

overdispersion=overdispersion^2;
VARCOV=0.5*(VARCOV+VARCOV'); %Impose symmetry 
if nargin ==6
if draw_type==0
    draw=posterior_mode;
elseif draw_type==1
    draw=mvnrnd(posterior_mode',VARCOV);
end
else
    draw=varargin{:}';
end
draws=zeros(number_of_draws,length(posterior_mode));
accepts=zeros(number_of_draws,1);
log_densities=zeros(number_of_draws,1);
candidates=zeros(number_of_draws,length(posterior_mode));


tic;
%%%%% Posterior Draws %%%%%
[log_density, lagged_expectation]=model_loglike(draw,selection_matrix,Data,T,K,model_name); 
prior=log_prior_function(draw); 
draws(1,:)=draw;
accepts(1,1)=1;
log_densities(1,1)=-log_density;
candidates(1,:)=draw;
lagged_expectations(1,1)=lagged_expectation;

mywaitbar = waitbar(0,'Running Metropolis Hastings');
for i=2:number_of_draws % Metropolis iterations
    accept=0;
    candidate = mvnrnd(draw',overdispersion*VARCOV);
    if min(candidate<=Pparameter(4,:))>0 && min(candidate>=Pparameter(3,:))>0
        prior_candidate = log_prior_function(candidate);
        [log_density_candidate, lagged_expectation_candidate]=model_loglike(candidate,selection_matrix,Data,T,K,model_name);
        r = exp(log_density-prior-log_density_candidate+prior_candidate);
        U=rand;
        if U<r;
            draw=candidate; accept=1; log_density=log_density_candidate; prior=prior_candidate; 
        end
        lagged_expectation=lagged_expectation_candidate;
    end
    draws(i,:)=draw;
    accepts(i,1)=accept;
    log_densities(i,1)=-log_density;
    candidates(i,:)=candidate;
    lagged_expectations(i,1)=lagged_expectation;
    proportion_finished=100*i/number_of_draws;
    current_acceptance_rate=mean(accepts(1:i,1));
    waitbar(proportion_finished/100,mywaitbar,{sprintf('Running Metropolis Hastings'); sprintf('Percent Finished: %2.2f. Current Acceptance Rate: %2.2f',proportion_finished,current_acceptance_rate)});
end
close(mywaitbar);
time_to_run=toc;
hours_to_run=floor(time_to_run/3600);
minutes_to_run=floor((time_to_run-hours_to_run*3600)/60);
seconds_to_run=time_to_run-hours_to_run*3600-minutes_to_run*60;


results = draws(burn_in+1:number_of_draws,:);
results_densities=log_densities(burn_in+1:number_of_draws,1);

display_results=num2str([mean(results)' median(results)' diag(cov(results)).^(1/2)]);
[junk width]=size(display_results);
display_results=[['Mean' blanks(floor((width-13)/2)) 'Median'  blanks(width-13-floor((width-13)/2))  'Std']; repmat('-',[1 width]); display_results; repmat('-',[1 width])];
disp(blanks(5)')
disp('Estimation Results:')
disp(display_results)
disp(['Final Acceptance Rate:' blanks(2) num2str(mean(accepts(burn_in+1:number_of_draws,1)))])
disp(['Execution Time:' blanks(2) num2str(hours_to_run) ' hour(s), ' num2str(minutes_to_run) ' minute(s), and '  num2str(seconds_to_run) ' second(s)'])
disp(['Execution Time per Iteration:' blanks(2) num2str(time_to_run/number_of_draws) ' second(s)'])

save(save_name,'results','draws','accepts','log_densities','candidates','VARCOV','posterior_mode','selection_matrix','Data','T','K','model_name','log_prior_function','Ptype','Pparameter')


