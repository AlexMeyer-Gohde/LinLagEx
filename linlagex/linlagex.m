%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%linlagex.m
%
%This file calls all other files needed to run the program
%
%One needs to have written an .m-file that calls this program to use the
%software (see examples)
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%15.01.2009 Fixed MA_Vector for spectral w/ no lagged exp. (fix is in
%spectral.m)
%
%added testing for significant imaginary elements in solution
%
%9.12.2009 Replaced the eval loop for systems without lagged expectations
%to set A_inf, etc. to zero matrices to save comp. time
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[num_eqs,num_endog]=size(A_0);%recovers number of endogenous variables
[num_eqs,num_exog]=size(G_0);%recovers number of exogenous variables
VARIABLE_NAMES=[ENDOGENOUS_VARIABLE_NAMES;EXOGENOUS_VARIABLE_NAMES];

if exist('it_max_value','var')==0
it_max_value=0;
end

if ischar(it_max_value)==1 %if the expectations reach into the infinite past (i.e. not an interger but infinity)
    limiting_matrix; %Runs the m file limiting_matrix to check whether the limiting matrices were declared and,...
                     %if not, to try and construct them symbolically
end

if ischar(it_max_value)==1 || it_max_value~=0 %if the expectations reach back into the infinite or finite past (i.e. there are past expectations)
matrix_function; %calls the m file of the same name to construct the anonymous functions for calculating the non-autonomous MA-recusions
end



if exist('growth_restriction','var')==0 %sets the growth restriction to 1 if not given by user
growth_restriction=1;
end
if exist('solution_method','var')==0 || strcmp(solution_method,'QZ')==0 && strcmp(solution_method,'AIM')==0 %If the solution method hasn't been chosen, or if a nonadmissible choice was made,
solution_method='QZ';%the default is QZ
end

if ischar(it_max_value)==1 %if the expectations reach into the infinite past (i.e. not an interger but infinity)
    if exist('tolerance','var')==0 %sets the tolerance if not given by user
        tolerance=1e-10;
    end
    max_it=solve_max_it(A_j,B_j,C_j,F_j,G_j,tolerance); %solves for the number of past expectations to be included before using the autonomous...
                                                        %recursion based on the limiting coefficients 
else %if the expectations reach into the finite past or are not present at all
    max_it=it_max_value; 
end

if max_it>0 %if there are past expectations
    if num_endog*num_eqs*max_it<500
    solve_ma; %calls the m file solve_ma to construct the block tridiagonal system and solve for the ma coefficients
    else
    solve_ma_alt;
    end
else    %if there are no past expectations
    A_inf=A_0;B_inf=B_0;C_inf=C_0;F_inf=F_0;G_inf=G_0;%with no past expectations, the limiting matrices (A_inf, etc.) don't change (A_inf=A_0, etc.)
    %calls the program QZ_solve or sp_solve to find the saddle-path recursion
    %Theta_i=ALPHA_ZS*Theta_(i-1)+BETA_ZS*(N^(i))
    if strcmp(solution_method,'QZ')==1
        [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=QZ_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,growth_restriction); %direct to QZ
    elseif strcmp(solution_method,'AIM')==1
        [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=sp_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,growth_restriction); %direct to AIM
    end
    MA_VECTOR=BETA_ZS;%needed for impulse.m
    LHS=[];%needed for anticipated impulse.m
    trisolve='';%needed for error checking. Since there is no tridiagonal system to solve, it is always successful
end

if isnumeric(existence_uniqueness)==1 && isempty(trisolve)==1 %If the asymptotic solution exists and is unique and if MATLAB solved the tri-diagonal system
    if exist('RUN_IMPULSE','var')==0 || RUN_IMPULSE==1 %If either the option has not been set or if the option RUN_IMPULSE has been set to 1
        impulse_response;%calculates impulse responses to innovation
    end
    if exist('RUN_ANTICIPATED_IMPULSE','var')==0 || RUN_ANTICIPATED_IMPULSE==1 %If either the option has not been set or if the option RUN_ANTICIPATED_IMPULSE has been set to 1
        anticipated_impulse;%calculates impulse responses to anticipated innovations
    end
    if exist('RUN_SPECTRAL','var')==0 || RUN_SPECTRAL==1 %If either the option has not been set or if the option RUN_SPECTRAL has been set to 1
        spectral;%calculates frequency domain moments
    end
    if exist('RUN_SIMULATION','var')==0 || RUN_SIMULATION==1 %If either the option has not been set or if the option RUN_SIMULATION has been set to 1
        simulation;%runs a simulation and calculates simulated moments
    end
%%%%%% Displays the applicable error message if the software was unable to
%%%%%% deliver a result
elseif strcmp(existence_uniqueness,'unstable')==1 
    disp('Sorry, your system is unstable with respect to the growth restriction')
    disp('')
elseif strcmp(existence_uniqueness,'indeterminate')==1
    disp('Sorry, your system is indeterminate with respect to the growth restriction')
    disp('')
elseif strcmp(existence_uniqueness,'non-translatable')==1
    disp('Sorry, I was unable to ''translate'' your eigenvalues')
    disp('')
elseif strcmp(existence_uniqueness,'imaginary')==1
    disp('Sorry, QZ/AIM led to significant imaginary entries in the solution')
    disp('')
elseif isempty(trisolve)==0
    disp('I was unable to solve the tri-diagonal system')
    disp('MATLAB returned the following error while trying to solve:')
    trisolve
end