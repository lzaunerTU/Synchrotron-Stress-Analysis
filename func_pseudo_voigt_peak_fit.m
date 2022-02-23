function [fitted_curve,gof] = func_pseudo_voigt_peak_function_fit(x,y,x0)
% Fits the x-y input data with a Pseudo-Voigt peak function (Originlab) and 
% returns the fitted parameters y0, xc, w, A, mu, as well as the goodness of 
% the fit. Start values for the fit parameters are input via x0.
% y0 ... vertical offset
% xc ... fitted peak position
% w  ... FWHM
% A  ... area under peak
% mu ... profile shape factor
% author:   Lukas Zauner
% contact:  lukas.zauner@tuwien.ac.at
% date:     Q4, 2021

    fitfun = fittype( @(y0,xc,w,A,mu,x) y0+A*[mu*2*w./(pi*4*(x-xc).^2+w^2)+(1-mu)*(sqrt(4*log(2))/(sqrt(pi)*w))*exp(-4*log(2).*((x-xc).^2)/w^2)]);
    [fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);
end

