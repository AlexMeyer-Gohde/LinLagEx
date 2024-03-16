%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%simulation.m
%
%This program borrows from
%simul.m VERSION 4.1, May 2003, COPYRIGHT H. UHLIG.
%for running a simulation and HP-filtering the simulation
%
%Moments are calculated (as in spectral.m) to allow for
%a broader (w.r.t. H. UHLIG's toolkit) range of cross-correlations 
%to be calculated
%
%THIS VERSION: 1.1 December 9, 2009
%Fixed 124, 126: 100 replaced with w sim length 9.1.2009
%Fixed 105: non HP-filtered sims demeaned 20.5.2009
%   Thanks go to Jan-Oliver Menz and Lena Vogel for pointing this out 
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma_1_sim=[];
Gamma_1_sim_hp=[];
Xcorr_sim=[];
Xcorr_sim_hp=[];

if exist('INITIAL_SIMULATION','var')==0
INITIAL_SIMULATION=100;
end
if exist('SIMULATION_LENGTH','var')==0
SIMULATION_LENGTH=100;
end
if exist('HP_LAMBDA','var')==0
if exist('PERIOD','var')==0
disp('You did not set PERIOD (number of periods per year)')
disp('I shall set HP_LAMBDA under the assumption that your model')
disp('is quarterly (four periods per year).')
PERIOD=4;
end
HP_LAMBDA=1600*(PERIOD/4)^4; %Ravn, Morten and Harald Uhlig (2001) "On Adjusting the HP-Filter for the Frequency of Observations"
end
if exist('CORRELATION_HORIZON','var')==0
CORRELATION_HORIZON=6;
end

SIMULATION_SHOCKS=chol(Omega,'lower')*randn(length(Omega),SIMULATION_LENGTH+INITIAL_SIMULATION);
SIM_PATH=zeros(num_endog,SIMULATION_LENGTH+INITIAL_SIMULATION);
SIM_PATH_EXOG=zeros(num_exog,SIMULATION_LENGTH+INITIAL_SIMULATION);

SIM_PATH_EXOG(:,1)=SIMULATION_SHOCKS(:,1);
for j=2:SIMULATION_LENGTH+INITIAL_SIMULATION
SIM_PATH_EXOG(:,j)=N*SIM_PATH_EXOG(:,j-1)+SIMULATION_SHOCKS(:,j);
end


if SIMULATION_LENGTH+INITIAL_SIMULATION<=max_it+1
    SIM=MA_VECTOR*SIMULATION_SHOCKS;
        for k=1:num_endog
            for j=1:SIMULATION_LENGTH+INITIAL_SIMULATION
                for i=1:j
                    SIM_PATH(k,j)=SIM_PATH(k,j)+SIM((j-i)*num_endog+k,i);
                end
            end
        end
elseif SIMULATION_LENGTH+INITIAL_SIMULATION>max_it+1
        Response=[MA_VECTOR];
    
        MA_VECTOR_LAST=MA_VECTOR(max_it*num_endog+1:(max_it+1)*num_endog,:);
    
    for j=max_it+2:SIMULATION_LENGTH+INITIAL_SIMULATION
        MA_VECTOR_LAST=ALPHA_ZS*MA_VECTOR_LAST+BETA_ZS*N^(j-1);
        Response=[Response;MA_VECTOR_LAST];
    end
    SIM=Response*SIMULATION_SHOCKS;
    for k=1:num_endog
        for j=1:SIMULATION_LENGTH+INITIAL_SIMULATION
            for i=1:j
                SIM_PATH(k,j)=SIM_PATH(k,j)+SIM((j-i)*num_endog+k,i);
            end
        end
    end
end
SIM_PATH=[SIM_PATH;SIM_PATH_EXOG];
SIM_PATH=SIM_PATH';

%%%The following is taken directly from simul.m VERSION 4.1, May 2003,
%%%COPYRIGHT H. UHLIG. and "is due to Gerard A. Pfann"
LENGTH = max(size(SIM_PATH));
   HP_mat = [1+HP_LAMBDA, -2*HP_LAMBDA, HP_LAMBDA,              sparse(1,LENGTH-3);
             -2*HP_LAMBDA,1+5*HP_LAMBDA,-4*HP_LAMBDA,HP_LAMBDA, sparse(1,LENGTH-4);
                           sparse(LENGTH-4,LENGTH);
              sparse(1,LENGTH-4),HP_LAMBDA,-4*HP_LAMBDA,1+5*HP_LAMBDA,-2*HP_LAMBDA;     
              sparse(1,LENGTH-3),          HP_LAMBDA,   -2*HP_LAMBDA, 1+HP_LAMBDA  ];
   for iiiii=3:LENGTH-2;
     HP_mat(iiiii,iiiii-2)=HP_LAMBDA;
     HP_mat(iiiii,iiiii-1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii)=1+6*HP_LAMBDA;
     HP_mat(iiiii,iiiii+1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii+2)=HP_LAMBDA;
   end;
SIM_PATH_tr=HP_mat\SIM_PATH;
SIM_PATH_hp=SIM_PATH-SIM_PATH_tr; 

SIM_PATH=SIM_PATH(INITIAL_SIMULATION+1:SIMULATION_LENGTH+INITIAL_SIMULATION,:);
SIM_PATH_MOMENTS=SIM_PATH-repmat(mean(SIM_PATH),[SIMULATION_LENGTH,1]); %%%%%5/20/2009
SIM_PATH_hp=SIM_PATH_hp(INITIAL_SIMULATION+1:SIMULATION_LENGTH+INITIAL_SIMULATION,:);

%%%%%%Begin moment calculations%%%%%%%%%%%
if exist('DISPLAY_STD','var')==0 || DISPLAY_STD==1
if max(abs(eig(N)))>=1 || max(abs(eig(ALPHA_ZS)))>=1
disp('simulation.m')
disp('Sorry, the system is not covariance stationary')
disp('I will only calculate HP-filtered moments')
NON_STATIONARY=1;
else
NON_STATIONARY=0;
Sample_Standard_Deviations=[char(VARIABLE_NAMES),repmat(char(32),[num_endog+num_exog,3]), num2str((diag((SIM_PATH_MOMENTS)'*(SIM_PATH_MOMENTS)/SIMULATION_LENGTH)).^(1/2),'% 0.5f')]
disp(' ')
end
Sample_Standard_Deviations_HP=[char(VARIABLE_NAMES),repmat(char(32),[num_endog+num_exog,3]), num2str((diag((SIM_PATH_hp)'*(SIM_PATH_hp)/SIMULATION_LENGTH)).^(1/2),'% 0.5f')]
disp(' ')
end

for j=0:CORRELATION_HORIZON
if NON_STATIONARY==0
Gamma_1_sim{j+1}=SIM_PATH_MOMENTS(1:SIMULATION_LENGTH-j,:)'*SIM_PATH_MOMENTS(1+j:SIMULATION_LENGTH,:);
end
Gamma_1_sim_hp{j+1}=SIM_PATH_hp(1:SIMULATION_LENGTH-j,:)'*SIM_PATH_hp(1+j:SIMULATION_LENGTH,:);
end
for j=1:2*length(Gamma_1_sim_hp)-1
for i=1:num_endog+num_exog
if j-length(Gamma_1_sim_hp)<0
if NON_STATIONARY==0
Xcorr_sim{i}(:,j)=Gamma_1_sim{abs(j-(length(Gamma_1_sim_hp)+1))}(:,i)./((Gamma_1_sim{1}(i,i).^(1/2))*(diag(Gamma_1_sim{1}).^(1/2)));
end
Xcorr_sim_hp{i}(:,j)=Gamma_1_sim_hp{abs(j-(length(Gamma_1_sim_hp)+1))}(:,i)./((Gamma_1_sim_hp{1}(i,i).^(1/2))*(diag(Gamma_1_sim_hp{1}).^(1/2)));
else
if NON_STATIONARY==0
Xcorr_sim{i}(:,j)=Gamma_1_sim{abs(j-(length(Gamma_1_sim_hp))+1)}(i,:)'./((Gamma_1_sim{1}(i,i).^(1/2))*(diag(Gamma_1_sim{1}).^(1/2)));
end
Xcorr_sim_hp{i}(:,j)=Gamma_1_sim_hp{abs(j-(length(Gamma_1_sim_hp))+1)}(i,:)'./((Gamma_1_sim_hp{1}(i,i).^(1/2))*(diag(Gamma_1_sim_hp{1}).^(1/2)));
end
end
end


%%%%%%Begin simulation plots%%%%%%%%%%%
if exist('PLOT_SIMULATION','var')==0 || PLOT_SIMULATION==1 %If either the option has not been set or if the option PLOT_SIMULATION has been set to 1
    if exist('SIMULATION_SELECT','var')==0
        SIMULATION_SELECT=1:num_endog+num_exog;
    end
    % If the user didn't specify the number of periods in a year, quarters are
    % assumed
    if exist('PERIOD','var')==0
        PERIOD=4;
    end

    TIME=(0:SIMULATION_LENGTH-1)/PERIOD;
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(SIMULATION_SELECT)-1);0.25, 0.25,0.25])
    SIM=plot(TIME,[SIM_PATH(:,SIMULATION_SELECT)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title({''Simulation'';''First %d Periods Discarded''})',INITIAL_SIMULATION))
    legend([VARIABLE_NAMES(SIMULATION_SELECT,1)],'Location','Best');
    ylabel('% Deviations from Steady-State');
    xlabel('Years');
    hold on
    plot(TIME, 0*TIME,'k')
    hold off
    pause;
    clf('reset')
    for i=1:length(VARIABLE_NAMES)
        subplot(ceil(length(VARIABLE_NAMES)^(1/2)),round(length(VARIABLE_NAMES)^(1/2)),i);
        SIM=plot(TIME,[SIM_PATH(:,i)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
        legend([VARIABLE_NAMES(i,1)],'Location','Best');
        hold on
        plot(TIME, 0*TIME,'k')
        hold off
        ylabel('% Deviations from Steady-State');
        xlabel('Years');
        eval(sprintf('title({''Simulation'';''First %d Periods Discarded''})',INITIAL_SIMULATION))
    end
    pause;
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(SIMULATION_SELECT)-1);0.25, 0.25,0.25])
    SIM=plot(TIME,[SIM_PATH_hp(:,SIMULATION_SELECT)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title({''Simulation (HP-Filtered)'';''First %d Periods Discarded''})',INITIAL_SIMULATION))
    legend([VARIABLE_NAMES(SIMULATION_SELECT,1)],'Location','Best');
    ylabel('% Deviations from Steady-State');
    xlabel('Years');
    hold on
    plot(TIME, 0*TIME,'k')
    hold off
    pause;
    clf('reset')
    for i=1:length(VARIABLE_NAMES)
        subplot(ceil(length(VARIABLE_NAMES)^(1/2)),round(length(VARIABLE_NAMES)^(1/2)),i);
        SIM=plot(TIME,[SIM_PATH_hp(:,i)]',':.','MarkerEdgeColor','auto','MarkerSize',8);
        legend([VARIABLE_NAMES(i,1)],'Location','Best');
        hold on
        plot(TIME, 0*TIME,'k')
        hold off
        ylabel('% Deviations from Steady-State');
        xlabel('Years');
        eval(sprintf('title({''Simulation (HP-Filtered)'';''First %d Periods Discarded''})',INITIAL_SIMULATION))
    end
    pause;
end

%%%%%%Begin correlation plots%%%%%%%%%%%
if exist('PLOT_CORRELATION','var')==0 || PLOT_CORRELATION==1 %If either the option has not been set or if the option PLOT_CORRELATION has been set to 1
if exist('CORRELATION_SELECT','var')==0
CORRELATION_SELECT=1:(num_endog+num_exog);
end
if NON_STATIONARY==0
for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_sim{CORRELATION_SELECT(j)}(CORRELATION_SELECT,1:2*CORRELATION_HORIZON+1)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''Sample Cross-Correlations of Selected Variables at t+j with %s at t'')',VARIABLE_NAMES{CORRELATION_SELECT(j)}))
    legend(VARIABLE_NAMES(CORRELATION_SELECT,1),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
end
for j=1:length(CORRELATION_SELECT)
    clf('reset')
    set(gcf,'DefaultAxesColorOrder',[hsv(length(CORRELATION_SELECT)-1);0.25, 0.25,0.25])
    CROSS=plot(-CORRELATION_HORIZON:CORRELATION_HORIZON,Xcorr_sim_hp{CORRELATION_SELECT(j)}(CORRELATION_SELECT,1:2*CORRELATION_HORIZON+1)',':.','MarkerEdgeColor','auto','MarkerSize',8);
    eval(sprintf('title(''HP-Filtered Sample Cross-Correlations of Selected Variables at t+j with %s at t'')',VARIABLE_NAMES{CORRELATION_SELECT(j)}))
    legend(VARIABLE_NAMES(CORRELATION_SELECT,1),'Location','Best');
    ylabel('Correlation Coefficient');
    xlabel('j');
    axis([-Inf,Inf,-1,1])
    pause;
end
end
