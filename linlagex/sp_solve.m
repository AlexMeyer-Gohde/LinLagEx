function [ALPHA_ZS,BETA_ZS,sorted_eigenvalues,existence_uniqueness]=...
    sp_solve(A_inf,B_inf,C_inf,F_inf,G_inf,N,num_eqs,num_endog,num_exog,uprbnd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%sp_solve.m
%
%
%The following program has been adapted from: 
%The Anderson-Moore Algorithm
%http://www.federalreserve.gov/pubs/oss/oss4/sp_solve.zip
%Version 6
%
%The following announcement conforms to the request by Gary Anderson
%regarding the AIM method
%---------------------------------------------------------
%Note: This program is in the public domain and may be used freely.
%however the authors would appreaciate acknowledgement of the source 
%by citation of any of the following papers
%
%Anderson, G and Moore, G
%"A Linear Algebraic Procedure for Solving Linear Perfect Foresight Models."
%Economics Letters, 17, 1985.
%
% Anderson, G.
%"Solving Linear Rational Expectations Models: A Horse Race"
%Computational Economics, 2008 vol 31 issue 2 pages 95-113
%
%Anderson, G.
%"A Reliable and Computationally Efficient Algorithm for Imposing the Saddle Point Property in Dynamic Models"
%Journal of Economic Dynamics and Control, Forthcoming
%---------------------------------------------------------
%
%This program uses the AIM algorithm to find the recursive solution for the
%autonomous recursion of MA-Coefficients.
%
%THIS VERSION: 1.1 December 9, 2009
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 nlag=1;
 nlead=1;
cof=[];
 scof=[];
 cofb=[];
 nlag=1;
 nlead=1;
 sorted_eigenvalues=[];
 existence_uniqueness=0;
 lgrts=[];
 ALPHA_ZS=[];
 BETA_ZS=[];



  % Numerical tolerances for aim


  condn  = 1.e-10;
  uprbnd = uprbnd + 1.e-6;

  % ---------------------------------------------------------------------
  % Construct structural coefficient matrix.
  % ---------------------------------------------------------------------



  % Construct cof matrix from cofg, cofh


  cof = [C_inf,zeros(num_endog,num_exog),B_inf,G_inf,A_inf,F_inf;zeros(num_exog,num_endog),-N,zeros(num_exog,num_endog),eye(num_exog,num_exog),zeros(num_exog,num_endog+num_exog)];

  [cofb,sorted_eigenvalues,ia,nex,nnum,lgrts,existence_uniqueness] = SPAmalg(cof,num_endog+num_exog,nlag,nlead,condn,uprbnd);

  if existence_uniqueness>1,
    disp(SPAimerr(existence_uniqueness));    
    if (existence_uniqueness==3) || (existence_uniqueness==35)  
        existence_uniqueness='unstable';
    elseif(existence_uniqueness==4) || (existence_uniqueness==45) 
        existence_uniqueness='indeterminate';
    else
        existence_uniqueness='non-translatable';
    end
  else
    ALPHA_ZS=cofb(1:num_endog,1:num_endog);
    psi=-G_inf-F_inf*N;
    phi=inv(B_inf+A_inf*ALPHA_ZS);
    temp=phi*psi;
    F=-phi*A_inf;
    vec_BETA_ZS=inv(eye(num_endog*num_exog)-kron(N',F))*temp(:);
    BETA_ZS=reshape(vec_BETA_ZS,num_endog,num_exog);
        if abs(imag([ALPHA_ZS BETA_ZS]))>1E-10
            existence_uniqueness='imaginary';
            ALPHA_ZS=[];
            BETA_ZS=[];
        else
            ALPHA_ZS=real(ALPHA_ZS);
            BETA_ZS=real(BETA_ZS);
        end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%The following programs have been taken from: 
%The Anderson-Moore Algorithm
%http://www.federalreserve.gov/pubs/oss/oss4/sp_solve.zip
%Version 6
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  [b,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
%                       SPAmalg(h,neq,nlag,nlead,condn,uprbnd)
%
%  Solve a linear perfect foresight model using the matlab eig
%  function to find the invariant subspace associated with the big
%  roots.  This procedure will fail if the companion matrix is
%  defective and does not have a linearly independent set of
%  eigenvectors associated with the big roots.
% 
%  Input arguments:
% 
%    h         Structural coefficient matrix (neq,neq*(nlag+1+nlead)).
%    neq       Number of equations.
%    nlag      Number of lags.
%    nlead     Number of leads.
%    condn     Zero tolerance used as a condition number test
%              by numeric_shift and reduced_form.
%    uprbnd    Inclusive upper bound for the modulus of roots
%              allowed in the reduced form.
% 
%  Output arguments:
% 
%    b         Reduced form coefficient matrix (neq,neq*nlag).
%    rts       Roots returned by eig.
%    ia        Dimension of companion matrix (number of non-trivial
%              elements in rts).
%    nexact    Number of exact shiftrights.
%    nnumeric  Number of numeric shiftrights.
%    lgroots   Number of roots greater in modulus than uprbnd.
%    aimcode     Return code: see function aimerr.

function [b,rts,ia,nexact,nnumeric,lgroots,aimcode] = ...
                        SPAmalg(h,neq,nlag,nlead,condn,uprbnd)

if(nlag<1 | nlead<1) 
    error('Aim_eig: model must have at least one lag and one lead.');
end

% Initialization.
nexact   = 0;
nnumeric = 0;
lgroots  = 0;
iq       = 0;
aimcode    = 0;
b=0;
qrows = neq*nlead;
qcols = neq*(nlag+nlead);
bcols = neq*nlag;
q        = zeros(qrows,qcols);
rts      = zeros(qcols,1);

% Compute the auxiliary initial conditions and store them in q.

[h,q,iq,nexact] = SPExact_shift(h,q,iq,qrows,qcols,neq);
   if (iq>qrows) 
      aimcode = 61;
      return;
   end

[h,q,iq,nnumeric] = SPNumeric_shift(h,q,iq,qrows,qcols,neq,condn);
   if (iq>qrows) 
      aimcode = 62;
      return;
   end

%  Build the companion matrix.  Compute the stability conditions, and
%  combine them with the auxiliary initial conditions in q.  

[a,ia,js] = SPBuild_a(h,qcols,neq);

if (ia ~= 0)
   [w,rts,lgroots] = SPEigensystem(a,uprbnd);
   q = SPCopy_w(q,w,js,iq,qrows);
end

test = nexact+nnumeric+lgroots;
   if (test > qrows) aimcode = 3;
   elseif (test < qrows) aimcode = 4;
   end

% If the right-hand block of q is invertible, compute the reduced form.

if(aimcode==0)
[nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,condn);


if ( nonsing & aimcode==0) aimcode =  1;
elseif (~nonsing & aimcode==0) aimcode =  5;
elseif (~nonsing & aimcode==3) aimcode = 35;
elseif (~nonsing & aimcode==4) aimcode = 45;
end
end;
%%
% [h,q,iq,nexact] = exact_shift(h,q,iq,qrows,qcols,neq)
%
% Compute the exact shiftrights and store them in q.

function [h,q,iq,nexact] = SPExact_shift(h,q,iq,qrows,qcols,neq)

hs=sparse(h);
nexact = 0;
left   = 1:qcols;
right  = qcols+1:qcols+neq;
zerorows = find( sum(abs( hs(:,right)' ))==0 );

while( any(zerorows) && iq <= qrows )
   nz = length(zerorows);
   q(iq+1:iq+nz,:) = hs(zerorows,left);
   hs(zerorows,:)   = SPShiftright(hs(zerorows,:),neq);
   iq     = iq + nz;
   nexact = nexact + nz;
   zerorows = find( sum(abs( hs(:,right)' ))==0 );
end
h=full(hs);
%%
%  [h,q,iq,nnumeric] = ...
%             SPNumeric_shift(h,q,iq,qrows,qcols,neq,condn)
%
% Compute the numeric shiftrights and store them in q.

function [h,q,iq,nnumeric] = SPNumeric_shift(h,q,iq,qrows,qcols,neq,condn)

nnumeric = 0;
left     = 1:qcols;
right    = qcols+1:qcols+neq;

[Q,R,E]  = qr( h(:,right) );
zerorows = find( abs(diag(R)) <= condn );

while( any(zerorows) && iq <= qrows )
   h=sparse(h);
   Q=sparse(Q);
   h = Q'*h;
   nz = length(zerorows);
   q(iq+1:iq+nz,:) = h(zerorows,left);
   h(zerorows,:)   = SPShiftright( h(zerorows,:), neq );
   iq       = iq + nz;
   nnumeric = nnumeric + nz;
   [Q,R,E] = qr( full(h(:,right)) );
   zerorows = find( abs(diag(R)) <= condn );
end
%%
function [y] = SPShiftright(x,n)

% [y] = shiftright(x,n)
%
%  Shift the rows of x to the right by n columns, leaving zeros in the
%  first n columns. 

[rows,cols] = size(x);

left  = 1:cols-n;
right = n+1:cols;

y = zeros(rows,cols);
y(:,right) = x(:,left);

return
%%
%  [a,ia,js] = SPBuild_a(h,qcols,neq)
%
%  Build the companion matrix, deleting inessential lags.
%  Solve for x_{t+nlead} in terms of x_{t+nlag},...,x_{t+nlead-1}.

function [a,ia,js] = SPBuild_a(h,qcols,neq)

left  = 1:qcols;
right = qcols+1:qcols+neq;
hs=sparse(h);
a0 = hs(:,right);
hs(:,left) = -hs(:,right)\hs(:,left);

%  Build the big transition matrix.

a = zeros(qcols,qcols);

if(qcols > neq)
   eyerows = 1:qcols-neq;
   eyecols = neq+1:qcols;
   a(eyerows,eyecols) = eye(qcols-neq);
end
hrows      = qcols-neq+1:qcols;
a(hrows,:) = hs(:,left);

%  Delete inessential lags and build index array js.  js indexes the
%  columns in the big transition matrix that correspond to the
%  essential lags in the model.  They are the columns of q that will
%  get the unstable left eigenvectors. 

js       = 1:qcols;
zerocols = sum(abs(a)) == 0;
while( any(zerocols) )
    a(:,zerocols) = [];
    a(zerocols,:) = [];
    js(zerocols)  = [];
    zerocols = sum(abs(a)) == 0;
  end
ia = length(js);
%%
%  [w,rts,lgroots] = SPEigensystem(a,uprbnd)
%
%  Compute the roots and the left eigenvectors of the companion
%  matrix, sort the roots from large-to-small, and sort the
%  eigenvectors conformably.  Map the eigenvectors into the real
%  domain. Count the roots bigger than uprbnd.

function [w,rts,lgroots] = SPEigensystem(a,uprbnd) 

[w,d]   = eig(a');
rts     = diag(d);
mag     = abs(rts);
[mag,k] = sort(-mag);
rts     = rts(k);

ws=sparse(w);
ws       = ws(:,k);

%  Given a complex conjugate pair of vectors W = [w1,w2], there is a
%  nonsingular matrix D such that W*D = real(W) + imag(W).  That is to
%  say, W and real(W)+imag(W) span the same subspace, which is all
%  that aim cares about. 

ws = real(ws) + imag(ws);

lgroots = sum(abs(rts) > uprbnd);

w=full(ws);
%%
% q = SPCopy_w(q,w,js,iq,qrows)
%
%  Copy the eigenvectors corresponding to the largest roots into the
%  remaining empty rows and columns js of q 

function  q = SPCopy_w(q,w,js,iq,qrows)

if(iq < qrows)
   lastrows = iq+1:qrows;
   wrows    = 1:length(lastrows);
   q(lastrows,js) = w(:,wrows)';
end
%%
% [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,b,condn);
%
% Compute reduced-form coefficient matrix, b.

function [nonsing,b] = SPReduced_form(q,qrows,qcols,bcols,neq,condn);
b=[];
qs=sparse(q);
left = 1:qcols-qrows;
right = qcols-qrows+1:qcols;
nonsing = rcond(full(qs(:,right))) > condn;
if(nonsing)
   qs(:,left) = -qs(:,right)\qs(:,left);
b = qs(1:neq,1:bcols);
b = full(b);
else  %rescale by dividing row by maximal qr element
disp('inverse condition number small, rescaling ')
themax=max(abs(qs(:,right)),[],2);
oneover = diag(1 ./ themax);
nonsing = rcond(full(oneover *qs(:,right))) > condn;
if(nonsing)
   qs(:,left) = -(oneover*qs(:,right))\(oneover*(qs(:,left)));
b = qs(1:neq,1:bcols);
b = full(b);
end
end
%%
function e = SPAimerr(c);
% e = aimerr(c);
%
% Interpret the return codes generated by the aim routines.

% The return code c = 2 is used by aim_schur.m but not by aim_eig.m.

    if(c==1)  e='Aim: unique solution.';
elseif(c==2)  e='Aim: roots not correctly computed by real_schur.';
elseif(c==3)  e='Aim: too many big roots.';
elseif(c==35) e='Aim: too many big roots, and q(:,right) is singular.';
elseif(c==4)  e='Aim: too few big roots.';
elseif(c==45) e='Aim: too few big roots, and q(:,right) is singular.';
elseif(c==5)  e='Aim: q(:,right) is singular.';
elseif(c==61) e='Aim: too many exact shiftrights.';
elseif(c==62) e='Aim: too many numeric shiftrights.';
else          e='Aimerr: return code not properly specified';
end

return