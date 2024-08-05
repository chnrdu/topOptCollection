function [ x, as, lmid ] = ocUpdate( loop, xT, dg0, g1, dg1, ocPar, xOld, xOld1, ...
    as, beta, restartAs )
% -------------------------------- definition of asymptotes and move limits
[xU,xL] = deal(min(xT+ocPar(1),1), max(xT-ocPar(1),0));
if (loop<2.5 || restartAs==1)
    as = xT+[-0.5,0.5].*(xU-xL)./(beta+1);
else
    tmp = (xT-xOld).*(xOld-xOld1);
    gm = ones(length(xT),1);
    [gm(tmp>0), gm(tmp<0)] = deal(ocPar(3),ocPar(2));
    as = xT+gm.*[-(xOld-as(:,1)),(as(:,2)-xOld)];
end
xL = max( 0.9*as(:,1)+0.1*xT,xL);                    % adaptive lower bound
xU = min( 0.9*as(:,2)+0.1*xT,xU);                    % adaptive upper bound
% ----- split (+) and (-) parts of the objective and constraint derivatives
p0_0 = (dg0>0).*dg0; q0_0 = (dg0<0).*dg0;
p1_0 = (dg1>0).*dg1; q1_0 = (dg1<0).*dg1;
[p0,q0] = deal(p0_0.*(as(:,2)-xT).^2,-q0_0.*(xT-as(:,1)).^2);
[p1,q1] = deal(p1_0.*(as(:,2)-xT).^2,-q1_0.*(xT-as(:,1)).^2);
% ---------------------------------------- define the primal projection map
primalProj = @(lm) min(xU,max(xL,(sqrt(p0+lm*p1).*as(:,1)+sqrt(q0+lm*q1).*as(:,2))...
    ./(sqrt(p0+lm*p1)+sqrt(q0+lm*q1))));
psiDual = @(lm) g1 - ( (as(:,2)-xT)'*p1_0 - (xT-as(:,1))'*q1_0 ) + ...
    sum(p1./(max(as(:,2)-primalProj(lm),1e-12)) + q1./(max(primalProj(lm)-as(:,1),1e-12)));
% ------------- compute the Lagrange multiplier and update design variables
lmUp = 1e6; x = xT; lmid = -1;
if psiDual( 0 ) * psiDual( lmUp ) < 0  % check if LM is within the interval
    lmid = fzero( psiDual, [ 0, lmUp ] );
    x = primalProj( lmid );
elseif psiDual( 0 ) < 0                       % constraint cannot be active
   lmid = 0;
   x = primalProj( lmid );
elseif psiDual( lmUp ) > 0                 % constraint cannot be fulfilled
   lmid = lmUp;
   x = primalProj( lmid );
end
%
end
%{
This code was written by F. Ferrari, J.K. Guest (a) and O. Sigmund (b)
(a) Johns Hopkins University, 21218 Baltimore MD (USA)
(b) Technical University of Denmark,2800 Lyngby (DK)
Please sent your comments to: fferrar3@jhu.edu

The code is intended for educational purposes and theoretical details
are discussed in the paper Ferrari,F. Sigmund,O., Guest,J.K. - Topology
Optimization with linearized buckling criteria in 250 lines of Matlab, SAMO
(2021). https://doi.org/10.1007/s00158-021-02854-x

The code as well as a postscript version of the paper can be
downloaded from: http://www.topopt.dtu.dk and www.ce.jhu.edu/topopt/

Disclaimer:
The authors reserves all rights but do not guarantee that the code is
free from errors. Furthermore, we shall not be liable in any event.
%}