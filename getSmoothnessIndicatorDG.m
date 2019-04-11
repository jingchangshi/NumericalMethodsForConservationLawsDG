function [se] = getSmoothnessIndicatorDG(u,iV,m)
% function [se] = getSmoothnessIndicatorDG(u,iV,m)
% Purpose: Compute nonlinear viscosity following Persson and Peraire (2006)

% Extract coefficients and compute smoothness measure
uh = iV*u;
S = uh(m+1,:).^2./sum(uh.*uh);
se = log(S);

end