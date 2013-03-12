function [lp alpha] = log_circ_vmpdf(alpha, thetahat, kappa)
%Adapted version of circ_vmpdf to give the log probability
%
% [p alpha] = circ_vmpdf(alpha, w, p)
%   Computes the circular von Mises pdf with preferred direction thetahat 
%   and concentration kappa at each of the angles in alpha
%
%   The vmpdf is given by f(phi) =
%   (1/(2pi*I0(kappa))*exp(kappa*cos(phi-thetahat)
%
%   Input:
%     alpha     angles to evaluate pdf at, if empty alphas are chosen to
%               100 uniformly spaced points around the circle
%     [thetahat preferred direction, default is 0]
%     [kappa    concentration parameter, default is 1]
%
%   Output:
%     p         von Mises pdf evaluated at alpha
%     alpha     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, FisherL
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu

% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen

%Modified to give log probability, Richard Mann (2010)
if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

%alpha = alpha(:);

if kappa < 1/0.04^2
% evaluate pdf
C = 1/(2*pi*besseli(0,kappa));
lp = log(C) + (kappa*cos(alpha-thetahat));
else %if kappa too large the above breaks, so use Gaussian limit
    lp = lognormpdf(alpha, thetahat, sqrt(1/kappa));
end
