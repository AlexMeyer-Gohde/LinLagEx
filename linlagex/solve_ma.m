%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%solve_ma.m
%
%Used when past expectations are present to construc and solve the block 
%tri-diagonal system
%
%This program uses the anonymous functions A_j,B_j,C_j,F_j and G_j and the
%matrices w/o lagged expectations (A_0,B_0,C_0,F_0 and G_0) to calculate
%matrix sums in order to form and solve the tri-diagonal system of linear 
%equations.
%
%The final restrictions are obtained through the QZ-Decomposition a la 
%Klein(2000) using the limiting matrix sums (given either by the user or
%solved for (1)(in the case of infinite lagged expectations) symbolically
%(done previously by the program limiting_matrix) or (2) (in the case of
%finite lagged expectations) by suming the matrices up to the highest
%lagged expectation included in the system.
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Predefine sparse system of equations
LHS=sparse((max_it+1)*num_eqs,(max_it+1)*num_endog);
RHS=zeros((max_it+1)*num_eqs,num_exog);
%Initial Entries
LHS(1:num_eqs,1:2*num_endog)=[B_0,A_0];
RHS(1:num_eqs,1:num_exog)=-(F_0*N+G_0);
%Begin Summation
A=A_0;
B=B_0;
C=C_0;
F=F_0;
G=G_0;



%Begin loop
for j=1:max_it-1

%Determine value of counter and sum matrices
A=A+A_j(j);
B=B+B_j(j);
C=C+C_j(j);
F=F+F_j(j);
G=G+G_j(j);
%enter these matrices into the block tridiagonal system
LHS(1+j*num_eqs:num_eqs+j*num_eqs,1+(j-1)*num_endog:(j+2)*num_endog)=[C,B,A];
RHS(1+j*num_eqs:num_eqs+j*num_eqs,1:num_exog)=-(F*(N^(j+1))+G*(N^j));
end

if ischar(it_max_value)==0 %If a finite number of lagged expectations are included
        for matj=[65,66,67,70,71]
                if eval(sprintf('eval([''exist(''''%s_inf'''')''])',char(matj)))==0%and if limiting matrix sums were not declared,
                    eval(sprintf('eval([''%s_inf=%s+%s_j(max_it);''])',char(matj),char(matj),char(matj)))%the sums are just the sums upto the
                                                                                                         %the last lagged epectation in the system
                end
        end
end
%calls the program QZ_solve or sp_solve to find the final set of restrictions
%Theta_max_it=ALPHA_ZS*Theta_(max_it-1)+BETA_ZS*(N^(max_it))
if strcmp(solution_method,'QZ')==1
        [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=QZ_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,growth_restriction); %direct to QZ
elseif strcmp(solution_method,'AIM')==1
        [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=sp_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,growth_restriction); %direct to AIM
end
if isnumeric(existence_uniqueness)==1
    %Enter Final restrictions into the block tridiagonal system

    LHS(1+max_it*num_eqs:(max_it+1)*num_eqs,1+(max_it-1)*num_endog:(max_it+1)*num_endog)=[-ALPHA_ZS,eye(num_endog,num_endog)];
    RHS(1+max_it*num_eqs:(max_it+1)*num_eqs,1:num_exog)=BETA_ZS*(N^(max_it));



    %Solve the sparse system for the first max_it coefficients

    lasterror('reset');
    lastwarn('');
    MA_VECTOR=LHS\RHS;
    MA_VECTOR=full(MA_VECTOR);
    trisolve=lasterror;
    trisolve=trisolve.message;
    if isempty(lastwarn)==0
        trisolve=lastwarn;
    end
end