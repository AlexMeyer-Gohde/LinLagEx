%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%impulse_response.m
%
%
%Calculates and plots impulse responses to shocks one percent
% in size to each one of the exogenous variables z(t), keeping all
% others at zero in the first period. 
%
% Response: the response of all variables x(t), y(t), z(t) to each shock.
%   Response is of size (m_states+n_endog+k_exog)*HORIZON.
%   Since Response is overwritten, each time a new shock is analyzed, 
%   the results are collected in 
% Resp_mat = [ Response to first shock
%              Response to second shock
%              ...                     ]"
%
%
% The basic structure of the program has been taken/adapted from: 
% VERSION 4.0, November 2002, COPYRIGHT H. UHLIG.
% IMPRESP.M
% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.
%
%Prof. Uhlig's program has been modified to fit the VARMA structure of the
% recursive solution
%
%Some, but by far not most, of the options available for Prof. Uhlig's
% toolkit have been integrated here
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If a horizon for impulse responses wasn't declared by the user, it is set
% to 40
if exist('HORIZON','var')==0
HORIZON=40;
end



%Sets initial empty matrices
Resp_mat = [];
Response= [];
Exog_MA_temp=[];
Exog_MA=[];

%calculates impulse responses of the exogenopus variables upto the largest
%of either HORIZON or max_it+1 (easier for subsequent calculations)
for j=1:max(HORIZON,max_it+1)
Exog_MA_temp=N^(j-1);
Exog_MA=[Exog_MA;Exog_MA_temp];
end


if HORIZON<=max_it+1 %If the impulse response falls wholly within the MA part of the solution
    for j=1:num_exog
        Response=reshape(MA_VECTOR(:,j),[num_endog,max_it+1]);
        Exog_response=reshape(Exog_MA(1:num_exog*(max_it+1),j),[num_exog,max_it+1]);
        Resp_mat=[Resp_mat;Response;Exog_response];
    end
    Resp_mat=Resp_mat(:,1:HORIZON);
elseif HORIZON>max_it+1 %Elseif the impulse response must be expanded beyond the MA part of the solution
        Response=[MA_VECTOR];
    
        MA_VECTOR_LAST=MA_VECTOR(max_it*num_endog+1:(max_it+1)*num_endog,:);
 
    for j=max_it+2:HORIZON
        MA_VECTOR_LAST=ALPHA_ZS*MA_VECTOR_LAST+BETA_ZS*N^(j-1);
        Response=[Response;MA_VECTOR_LAST];
    end
    for j=1:num_exog
        Response_temp=reshape(Response(:,j),[num_endog,HORIZON]);
        Exog_response=reshape(Exog_MA(:,j),[num_exog,HORIZON]);
        Resp_mat=[Resp_mat;Response_temp;Exog_response];
    end
end

%%%%%%%%%%%%%%%%%Plots for impulse responses%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('PLOT_IMPULSE','var')==0 || PLOT_IMPULSE==1 %If either the option has not been set or if the option PLOT_IMPULSE has been set to 1
% If the user didn't specify a subset of variables to be displayed, all
% variables are displayed
if exist('IMPULSE_SELECT','var')==0
IMPULSE_SELECT=1:(num_endog+num_exog);
end

% If the user didn't specify the number of periods in a year, quarters are
% assumed
if exist('PERIOD','var')==0
PERIOD=4;
end

TIME=(-1*PERIOD:HORIZON-1)/PERIOD;
for j=1:num_exog
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(IMPULSE_SELECT)-1);0.25, 0.25,0.25])
    IMPULSE=plot(TIME,[zeros(length(IMPULSE_SELECT),PERIOD) Resp_mat((j-1)*(num_endog+num_exog)+IMPULSE_SELECT,:)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''Impulse Responses to a Shock in %s'')',EXOGENOUS_VARIABLE_NAMES{j}))
    legend(VARIABLE_NAMES(IMPULSE_SELECT,1),'Location','Best');
    ylabel('% Deviations from Steady-State');
    xlabel('Years since Shock Realization');
    hold on
    plot(TIME, 0*TIME,'k')
    hold off
    pause;
    clf('reset')
    for i=1:length(VARIABLE_NAMES)
        subplot(ceil(length(VARIABLE_NAMES)^(1/2)),round(length(VARIABLE_NAMES)^(1/2)),i); plot(TIME, 0*TIME,'k', TIME, [zeros(1,PERIOD) Resp_mat((j-1)*(num_endog+num_exog)+i,:)],'k:.','MarkerEdgeColor','auto','MarkerSize',8);
        legend(VARIABLE_NAMES(i,1),'Location','Best');
        ylabel('% Deviations from Steady-State');
        xlabel('Years since Shock Realization');
        eval(sprintf('title({''Impulse Response to a Shock in'';''%s''})',EXOGENOUS_VARIABLE_NAMES{j}))
    end
    pause;
end
end