function [detW, detV, R_p]= convergence(mixed_draws, number_of_chains,burn_in,skip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%convergence.m
%
%This program:
%               calculates and plots some multivariate convergence
%               diagnostics from Brooks, S.P and Gelman, A. (1998) "General
%               Methods for Monitoring Convergence of Iterative
%               Simulations,' Journal of Computational and Graphical
%               Statistics, Vol. 7, No. 4, Decc. 1998, pp. 434-455
%
%Inputs:        mixed_draws: matrix of stacked draws
%               number_of_chains: number of chains
%               burn_in: how many initial draws to discard
%               skip: how many draws to skip over. I.e if skip =10, then
%               the final matrix of draws to use for the convergence test=
%               row(1), row(11), row(21), ...
%
%Output:        R_p: a vector of Multivariate Scale Reduction factors as a
%                   function of the iterations (Brooks and Gelman p. 446).
%                   Should be close to 1 if convergence
%               detW: determinant of within variance. Should stabilize as a
%                   function of the number of iteration
%               detV: determinant of pooled variance. Should stabilize as a
%                   function of the number of iteration
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


[n p]=size(mixed_draws);
n=n/number_of_chains;
select=1:skip:n;
total_select=[];
for j=1:number_of_chains
total_select=[total_select select+(j-1)*n];
end
mixed_draws=mixed_draws(total_select,:);
[n p]=size(mixed_draws);
n=n/number_of_chains;

for i=2:n
Wtotal=zeros(p,p);
for j=1:number_of_chains
phi_j_dot(j,:)=mean(mixed_draws((j-1)*n+1:(j-1)*n+i,:));
W=cov(mixed_draws((j-1)*n+1:(j-1)*n+i,:));
Wtotal=Wtotal+W;
end
Wtotal=Wtotal/number_of_chains;
Bntotal=cov(phi_j_dot);
V_hat=(i-1)*Wtotal/i+(1+1/number_of_chains)*Bntotal;
detW(i)=det(Wtotal);
detV(i)=det(V_hat);
if cond(Wtotal)<1e10
R_p(i)=(i-1)/i+(number_of_chains+1)*max(eig(Wtotal\Bntotal))/number_of_chains;
else
R_p(i)=Inf;
end
end
R_p=R_p(burn_in+1:end);
detW=detW(burn_in+1:end);
detV=detV(burn_in+1:end);
figure
plot(1:n-burn_in, detW, ':k', 1:n-burn_in, detV, '-k')
legend('Determinant of Within Variance','Determinant of Pooled Variance')
figure
plot(1:n-burn_in, R_p, '-k')
legend('Multivariate Potential Scale Reduction Factor')