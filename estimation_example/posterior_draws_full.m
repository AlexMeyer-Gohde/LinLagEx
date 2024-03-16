%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%posterior_draws_full.m
%
%This program:
%               is a script to get and assess draws from the posterior
%               distribution using five independent chains for the
%               non-truncated model
%
%      NOTE: The file posterior_mode_full.mat, containing the
%      results from a numerical maximization of the log posterior
%      likelihood function, must already exist. If this is not the case, run
%      the file posterior_mode_full first.
%
%       Simply change to the directory where this file is located, type
%       posterior_draws_full, wait a couple of hours, and the results will be 
%       displayed and then saved in the file
%       mixed_posteriors_full.mat
%
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
number_of_chains=5;
number_of_draws=50000;
burn_in=30000;
overdispersion=0.8;
skip=10;
credible_set_percentile=0.05;
mode_estimate_data='posterior_mode_full.mat';
save_final_name='mixed_posteriors_full';
parameter_names={'$\lambda$'
'$\xi$'
'$\frac{1}{a_1}$'
'$\phi_{\pi}$'
'$\phi_y$'
'$\rho_{is}$'
'$\rho_{pc}$'
'$\phi_{R}$'
'$\sigma^{eps}_{is}$'
'$\sigma^{eps}_{pc}$'
'$\sigma^{eps}_{mp}$'};
bayesian_plots_and_correlations_model_name='input_example_me_bayesian_correlations';
CORRELATION_HORIZON=10;
IMPULSE_HORIZON=40;
variable_selection_matrix=[1 2 3];
ENDOGENOUS_VARIABLE_NAMES={'Inflation'
'Output Gap'
'Nominal Interest Rate'};
EXOGENOUS_VARIABLE_NAMES={'IS Shock'
'PC Shock'
'MP Shock'};
PERIOD=4;
use_this_many_posteriors=5000;
correlation_credible_set_percentiles=0.1;

decomp_credible_set_percentiles=0.1;
decomp_horizon=[1 2 5 10 15];


addpath('..\estimation_engine')
addpath('..\estimation_engine\sims_opt')
addpath('..\linlagex')

load(mode_estimate_data)

draw_type=0;
save_name='posterior_full_draw_1';
initial_draw=posterior_mode;
[results,draws,results_densities,lagged_expectations]=metropolis_hastings(mode_estimate_data, draw_type, number_of_draws, burn_in, overdispersion,save_name);
mixed_results=results;
mixed_draws=draws;
mixed_densities=results_densities;
mixed_lagged_expectations=lagged_expectations;

if number_of_chains>1
draw_type=1;
for j=2:number_of_chains
eval(sprintf('eval([''save_name=''''posterior_full_draw_%s'''';''])', num2str(j)))
[results,draws,results_densities,lagged_expectations]=metropolis_hastings(mode_estimate_data, draw_type, number_of_draws, burn_in, overdispersion,save_name);
mixed_results=[mixed_results; results];
mixed_draws=[mixed_draws; draws];
mixed_densities=[mixed_densities; results_densities];
mixed_lagged_expectations=[mixed_lagged_expectations;lagged_expectations];
end
[detW, detV_hat, R_p]=convergence(mixed_draws,number_of_chains,burn_in/skip,skip);
end
plots(mixed_results,Pparameter,Ptype,parameter_names)

[kept_draws number_of_parameters]=size(mixed_results);
for j=1:number_of_parameters
temp= sort(mixed_results(:,j));
credible_set(j,1)=temp(round(credible_set_percentile*kept_draws+1));
credible_set(j,2)=temp(round((1-credible_set_percentile)*kept_draws));
end

mean_final_results=mean(mixed_results);
median_final_results=median(mixed_results);
std_final_results=var(mixed_results).^(1/2);
marginal = marginal_density(mixed_results,mixed_densities);

display_results=num2str([mean_final_results' median_final_results' std_final_results' credible_set]);
[junk width]=size(display_results);
display_results=[['Mean' blanks(floor((width-20)/4)) 'Median' blanks(floor((width-20)/4)) 'Std' blanks(floor((width-20)/4)) num2str(100*credible_set_percentile) ' %'  blanks(width-20-3*floor((width-20)/4)) num2str(100*(1-credible_set_percentile)) ' %']; repmat('-',[1 width]); display_results; repmat('-',[1 width])];
disp(blanks(5)')
disp('Estimation results:')
disp(display_results)
disp(['Number of Chains:' blanks(2) num2str(number_of_chains)])
disp(['Total Draws per Chain:' blanks(2) num2str(number_of_draws)])
disp(['Burn In:' blanks(2) num2str(burn_in)])
disp(['Total Number of Mixed Draws:' blanks(2) num2str(number_of_chains*(number_of_draws-burn_in))])
disp(['Using Mode Estimate Data:' blanks(2) mode_estimate_data])
disp(['Log Marginal Density:' blanks(2) num2str(marginal)])

disp(blanks(5)')
disp(['Calculating Posterior Correlations and Impulses... Please wait.'])
[mean_correlations mean_impulses credible_corrs credible_imps median_corrs median_imps data_correlations]=bayesian_impulses_and_correlations(mixed_results, reshape(Data, [length(variable_selection_matrix) length(Data)/length(variable_selection_matrix)])', bayesian_plots_and_correlations_model_name, variable_selection_matrix,ENDOGENOUS_VARIABLE_NAMES, EXOGENOUS_VARIABLE_NAMES, CORRELATION_HORIZON, IMPULSE_HORIZON, PERIOD, use_this_many_posteriors,correlation_credible_set_percentiles);
disp(['Calculating Posterior Variance Decompositions... Please wait.'])
[mean_decomps median_decomps credible_decomps]=bayesian_variance_decomposition(mixed_results, bayesian_plots_and_correlations_model_name, variable_selection_matrix, ENDOGENOUS_VARIABLE_NAMES, EXOGENOUS_VARIABLE_NAMES, decomp_horizon, PERIOD, use_this_many_posteriors,decomp_credible_set_percentiles);
disp(blanks(2)')
disp('Mean Variance Decomposition:')
disp(num2str(cell2mat(mean_decomps')))
disp(blanks(2)')
disp('Credible Variance Decomposition:')
disp(num2str(cell2mat(credible_decomps')))

disp(['Results Saved In File:' blanks(2) mode_estimate_data])
save(save_final_name,'marginal','median_decomps','mean_decomps','credible_decomps','median_corrs','median_imps','mean_correlations','mean_impulses','credible_corrs','credible_imps', 'data_correlations','mixed_lagged_expectations','mixed_densities','mixed_results','mixed_draws','mean_final_results', 'median_final_results', 'std_final_results', 'VARCOV','posterior_mode','selection_matrix','Data','T','K','model_name','log_prior_function','Ptype','Pparameter')