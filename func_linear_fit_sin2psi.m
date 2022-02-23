function [k,intersect,rsquare] = func_linear_fit_sin2psi(sin2psi,d,NaNidx,nrofscans,dmin)
% Fits the sin2psi - d input data with a linear function and returns the
% fitted parameters k (slope), intersect (intersection of linear
% function with y axis), and the goodness of the fit.
% author:   Lukas Zauner
% contact:  lukas.zauner@tuwien.ac.at
% date:     Q4, 2021

    k = zeros(nrofscans,1);
    intersect = zeros(nrofscans,1);
    rsquare = zeros(nrofscans,1);
    
    for i = 1:nrofscans
        cur_sin2psi = sin2psi(NaNidx(:,i));                                 % get psi values where corresponding lattice spacing is not NaN because of a bad peak fitting result (low R2 value)
        cur_d = d(:,i);                                                     % get all d values in column i
        cur_d = cur_d(NaNidx(:,i));                                         % exclude all NaN values from current d values vector
        
        if length(cur_d) > dmin                                             % minimum number of d spacing necessary to allow for linear fitting
            [cur_fit,cur_gof] = fit(cur_sin2psi,cur_d,'poly1');             % linear fit to current d vs. sin2psi data
            k(i) = cur_fit.p1;                                              % get the slope of the current fit
            intersect(i) = cur_fit.p2;                                      % get the intersection of current fit with y-axis
            rsquare(i) = cur_gof.rsquare;                                   % get R2 value of current fit
        else
            k(i) = NaN;                                                     % if less than dmin values are in current d vector return NaN for all fitting parameters
            intersect(i) = NaN;
            rsquare(i) = NaN;
        end    
    end    
end

