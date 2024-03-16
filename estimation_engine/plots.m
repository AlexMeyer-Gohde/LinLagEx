function []=plots(results,Pparameter,Ptype,parameter_names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%plots.m
%
%This program:
%               produces histogram plots of parameter draws, as well as
%                   density plots for priors and posteriors
%                   -For posteriors, the densities are calculated using
%                   ksdensity, a kernel smoothing density estimator, using
%                   a normal kernel
%
%Inputs:        results : draws from metropolis
%               Pparameter: matrix containing, in this order,
%                the row vector of means, the row vector of standard deviations, the row
%                vector of maximal values and the row vector of minimal values
%               Ptype: cell array of prior distribution types defined as
%                strings (see create_log_priors.m for more information)
%               parameter_names: a cell array with strings containing
%                   parameter names
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[junk parameters]=size(results);
figure;
for i=1:parameters
    subplot(parameters,1,i), hist(results(:,i),50);
    lg=legend(parameter_names{i});
    set(lg, 'Interpreter','latex')
end


for i=1:parameters
[p_f,p_x]= ksdensity(results(:,i));
min_B=min(results(:,i));
max_B=max(results(:,i));
x=min_B:(max_B-min_B)/100:max_B;
    if strcmp(Ptype{i}, 'norm')
    type='normpdf';
    A=Pparameter(1,i);
    B=Pparameter(2,i);
    eval(sprintf('eval([''y=%s(x,%.20e,%.20e);''])',type,A,B));
    end
    if strcmp(Ptype{i}, 'beta')
    temp1=(Pparameter(1,i)-Pparameter(3,i))/(Pparameter(4,i)-Pparameter(3,i));
    temp2=Pparameter(2,i)/(Pparameter(4,i)-Pparameter(3,i));
    A=((1-temp1)*temp1^2)/(temp2^2)-temp1;
    B=A*(1/temp1-1);
    upper=Pparameter(4,i);
    lower=Pparameter(3,i);
    scalar=1/(beta(A,B)*((upper-lower)^(A+B-1)));
    eval(sprintf('eval([''y=%.20e*((x-%.20e).^(%.20e-1).*(%.20e-x).^(%.20e-1));''])',scalar,lower,A,upper,B));
    end
    if strcmp(Ptype{i}, 'gamma')
    type='gampdf';
    A=((Pparameter(1,i)-Pparameter(3,i))/Pparameter(2,i))^2;
    B=(Pparameter(1,i)-Pparameter(3,i))/A;
    lower=Pparameter(3,i);
    eval(sprintf('eval([''y=%s(x-%.20e,%.20e,%.20e);''])',type,lower,A,B));
    end
    if strcmp(Ptype{i}, 'invgamma1')
    nu = sqrt(2*(2+Pparameter(1,i)^2/Pparameter(2,i)^2));
    nu2 = 2*nu;
    nu1 = 2;
    err = 2*Pparameter(1,i)^2*gamma(nu/2)^2-(Pparameter(2,i)^2+Pparameter(1,i)^2)*(nu-2)*gamma((nu-1)/2)^2;
    while abs(nu2-nu1) > 1e-12
        if err > 0
            nu1 = nu;
            if nu < nu2
                nu = nu2;
            else
                nu = 2*nu;
                nu2 = nu;
            end
        else
            nu2 = nu;
        end
        nu =  (nu1+nu2)/2;
        err = 2*Pparameter(1,i)^2*gamma(nu/2)^2-(Pparameter(2,i)^2+Pparameter(1,i)^2)*(nu-2)*gamma((nu-1)/2)^2;
    end
    s = (Pparameter(2,i)^2+Pparameter(1,i)^2)*(nu-2);
    scalar=2/(gamma(nu/2)*(2/s)^(nu/2));
    eval(sprintf('eval([''y=%.20e./(x.^(%.20e+1).*exp(%.20e./(2*x.^2)));''])',scalar,nu,s));
    end
    if strcmp(Ptype{i}, 'invgamma2')
    A=2+(Pparameter(1,i)^2)/(Pparameter(2,i)^2);%%% 14.9.2009
    B=Pparameter(1,i)*(A-1);
    scalar=(B^A)/gamma(A);
    prior{1,i}=sprintf('-(1+%.20e)*log(x(%d))-%.20e/x(%d)',A,i,B,i);
    eval(sprintf('eval([''y=%.20e./(x.^(%.20e+1).*exp(%.20e./x));''])',scalar,A,B));
    end
    if strcmp(Ptype{i}, 'uniform')
    type='unifpdf';
    A=Pparameter(3,i);
    B=Pparameter(4,i);
    eval(sprintf('eval([''y=%s(x,%.20e,%.20e);''])',type,A,B));
    end
figure;
hold on
plot(p_x,p_f,'r','linewidth',2)
plot(x,y,'k','linewidth',2)
legend('Posterior', 'Prior')
title(['Prior and Posterior of ' parameter_names{i}],'Interpreter','LaTex')
hold off
end

