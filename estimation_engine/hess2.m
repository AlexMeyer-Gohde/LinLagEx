

% function hessian_mat = hessian(func,x,varargin)
% Computes second order partial derivatives
%
% INPUTS
%    func:           name of the function
%    x:              vector of variables around which the Hessian is calculated
%    varargin:       list of arguments following x
%
% OUTPUTS
%    hessian_matrix: Hessian matrix
%
% ALGORITHM
% Uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
%
% SPECIAL REQUIREMENTS
%    none
%  

% Copyright (C) 2001-2007 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

gstep=1e-2;
  %func = str2func(func);
  n=size(x,1);
  %h1=max(abs(x),options_.gstep*ones(n,1))*eps^(1/3);
  h1=max(abs(x),sqrt(gstep)*ones(n,1))*eps^(1/6);
  h_1=h1;
  xh1=x+h1;
  h1=xh1-x;
  xh1=x-h_1;
  h_1=x-xh1;
  xh1=x;
  f0=objective_function(x);
  f1=zeros(size(f0,1),n);
  f_1=f1;
  for i=1:n    
    xh1(i)=x(i)+h1(i);
    f1(:,i)=objective_function(xh1);
    xh1(i)=x(i)-h_1(i);
    f_1(:,i)=objective_function(xh1);
    xh1(i)=x(i);
    i=i+1;
  end
  xh_1=xh1;
  hessian_mat = zeros(size(f0,1),n*n);
  for i=1:n    
    if i > 1        
      k=[i:n:n*(i-1)];
      hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
    end     
    hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
    temp=f1+f_1-f0*ones(1,n);
    for j=i+1:n        
      xh1(i)=x(i)+h1(i);
      xh1(j)=x(j)+h_1(j);
      xh_1(i)=x(i)-h1(i);
      xh_1(j)=x(j)-h_1(j);
      hessian_mat(:,(i-1)*n+j)=-(-objective_function(xh1)-objective_function(xh_1)+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
      xh1(i)=x(i);
      xh1(j)=x(j);
      xh_1(i)=x(i);
      xh_1(j)=x(j);
      j=j+1;
    end    
    i=i+1;
  end 
  % 11/25/03 SA Created from Hessian_sparse (removed sparse)