function [log_prior_function]=create_log_priors(Pparameter, Ptype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%create_log_priors.m
%This program:
%               -creates a log prior density function for the vector of
%               parameters to be estimated
%
%Input:         Pparameter: matrix containing, in this order,
%                the row vector of means, the row vector of standard deviations, the row
%                vector of maximal values and the row vector of minimal values
%               Ptype: cell array of prior distribution types defined as
%                strings
%               -the following distributions are included:
%               (1) Normal (Ptype='norm')
%               (2) Beta (Ptype='beta')
%               (3) Gamma (Ptype='gamma')
%               (4) Inverse Gamma (Ptype='invgamma1' or 'invgamma2')
%                      invgamma1: % X ~ IG1(s,nu)
%                       X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
%
%                      invgamma2:X ~ IG2(s,nu)
%                       X = inv(Z) where Z ~ G(nu/2,2/s) (Gamma distribution)  
%               See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more
%               details.
%                   -Note that the hyper parameters are formed from the
%                   means, standard deviations, and bounds and not directly inputed
%                   by the user
%               (5) Uniform (Ptype='uniform' note that the upper and lower
%                       bounds are simply the largest and smallest of the integers
%                       entered into Pmean and Pvar when using the uniform
%                       distribution)
%
%Output:        log_prior_function: a function that returns the log prior density
%                   for a given parameter vector 
%
%THIS VERSION: 1.1 December 9, 2009
%
%
%
%  This file draws heavily on the .m files: draw_prior_density, inverse_gamma_specification, lpdfbeta, lpdfgam,
%  lpdfgbeta, lpdfig1, lpdfig2, lpdfnorm from dynare_v3 and dynare_4_0_4.
%
%  The primary difference is that this file creates a function on the desktop with all the
%  constants already calculated, circumventing the repetitious calculations
%  in the dynare source code - the savings in computational time is
%  essentially irrelevant, as the likelihood calculations are much more computationally intensive.
%
% Accordingly:
%   Copyright (C) 2003-2008 Dynare Team
%
%   This file is part of Dynare.
%
%   Dynare is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   Dynare is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
% A copy of the GNU General Public License associated with Dynare can be
% found in the directory /dynare_gpl_fdl.  
% If not, see <http://www.gnu.org/licenses/>.
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Pparameter,2)~=length(Ptype)
disp('Pparameter and Ptype are of different dimensions!')
disp('Make sure you have a mean, standard deviation, bounds and type for each parameter!')
else
scalar=0;
for i=1:size(Pparameter,2)
if strcmp(Ptype{i}, 'norm')
A=Pparameter(1,i);
B=Pparameter(2,i);
scalar=scalar-log(B)-log(2*pi)/2-log(normcdf(Pparameter(4,i))-normcdf(Pparameter(3,i)));
temp=2*(B)^2;
prior{1,i}=sprintf('-((x(%d)-%.20e)^2)/%.20e',i,A,temp);
end
if strcmp(Ptype{i}, 'beta')
temp1=(Pparameter(1,i)-Pparameter(3,i))/(Pparameter(4,i)-Pparameter(3,i));
temp2=Pparameter(2,i)/(Pparameter(4,i)-Pparameter(3,i));
A=((1-temp1)*temp1^2)/(temp2^2)-temp1;
B=A*(1/temp1-1);
upper=Pparameter(4,i);
lower=Pparameter(3,i);
scalar=scalar-betaln(A,B)-(A+B-1)*log(upper-lower);
prior{1,i}=sprintf('+(%.20e-1)*log(x(%d)-%.20e)+(%.20e-1)*log(%.20e-x(%d))',A,i,lower,B,upper,i);
end
if strcmp(Ptype{i}, 'gamma')
A=((Pparameter(1,i)-Pparameter(3,i))/Pparameter(2,i))^2;
B=(Pparameter(1,i)-Pparameter(3,i))/A;
lower=Pparameter(3,i);
scalar=scalar-gammaln(A)-A*log(B)+lower/B;
prior{1,i}=sprintf('+(%.20e-1)*log(x(%d)-%.20e)-x(%d)/%.20e',A,i,lower,i,B);
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
scalar=scalar+log(2)-gammaln(nu/2)-(nu/2)*log(2/s);
prior{1,i}=sprintf('-(%.20e+1)*log(x(%d))-%.20e/(2*(x(%d)^2))',nu,i,s,i);
end
if strcmp(Ptype{i}, 'invgamma2')
A=2+(Pparameter(1,i)^2)/(Pparameter(2,i)^2);%%% 14.9.2009
B=Pparameter(1,i)*(A-1);
scalar=scalar-gammaln(A)+A*log(B);
prior{1,i}=sprintf('-(1+%.20e)*log(x(%d))-%.20e/x(%d)',A,i,B,i);
end
if strcmp(Ptype{i}, 'uniform')
A=Pparameter(3,i);
B=Pparameter(4,i);
scalar=scalar-log(B-A);
prior{1,i}=sprintf('');
end
end
end
eval(sprintf('eval([''log_prior_function=@(x)%s;''])',strcat(num2str(scalar),prior{:})))

