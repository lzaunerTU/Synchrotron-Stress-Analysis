% author:   Lukas Zauner
% contact:  lukas.zauner@tuwien.ac.at, Institute of Materials Science and
%           Technology, TU Wien, Austria
% date:     Q4, 2021

%% This script extracts the stress state from synchrotron diffraction patterns (Debye-Scherrer rings) using the sin2(psi)-method
% The script uses 10 .dat files as input which result from azimuthal 
% integration of the ring pattern in 10° cakes between psi = 0 and 90° 
% using the DPDAK open-source software [1,2]. The script can handle multiple 
% scans at once (i.e., n by m array). If more than one scan should be
% evaluated, all corresponding results from azimuthal integration need to 
% be in the same .dat file, starting with the upper left corner and going 
% column wise from left to right (i.e., 4 x 4 array of scans --> every .dat
% file contains 32 columns [16 identical columns with radial diffraction 
% angle data, 16 columns with the intensity data from cake integration]).

% [1] G. Benecke, W. Wagermaier, C. Li, M. Schwartzkopf, G. Flucke, 
%       R. Hoerth, I. Zizak, M. Burghammer, E. Metwalli, P. Müller-Buschbaum, 
%       M. Trebbin, S. Förster, O. Paris, S. V. Roth, and P. Fratzl, 
%       J. Appl. Crystallogr. 47, 1797 (2014)

% [2] J. Keckes, M. Bartosik, R. Daniel, C. Mitterer, G. Maier, W. Ecker, 
%       J. Vila-Comamala, C. David, S. Schoeder, and M. Burghammer, 
%       Scr. Mater. 67, 748 (2012)

close all;
clear all;

%% Input of basic parameters
material = 'xxx';                                                           % material, just for reference
dir_data = 'C:\Users\...\';                                                  
                                                                            % directory to folder with .dat files (files need to be ordered based on increasing azimuthal angle)

lambda_sync = 0.82656;                                                      % [A] synchrotron X-ray wavelength
 
scan_columns = 2;                                                           % numbers of columns for evaluation (starting left)
scan_rows = 20;                                                             % numbers of rows for evaluation (starting on top)
nrofscans = scan_columns*scan_rows;                                         % number of scans for evaluation (can be smaller that total number of scans because only specified rows & columns are evaluated)
nrofscans_total = 400;                                                      % total number of scans processed in DPDAK and included in the .dat file
scan_stepsize = 250;                                                        % [nm] distance between scans, assumed to be equal in all directions

half_s2 = 2.575*10^(-3);                                                    % [GPa^-1] X-ray elastic constant 1/2s2 for material: see Birkholz - Thin Film Analysis by X-Ray Scattering (2006), p.259
psi_stress_free = 42;                                                       % stress free angle psi* for material, see Birkholz - Thin Film Analysis by X-Ray Scattering (2006), p.259


%% Input of peak fitting parameters
y0 = 0.1;                                                                   % background level or y offset
xc = 22.7;                                                                  % peak center
w = 0.1;                                                                    % FWHM
A = 70;                                                                     % area under peak
mu = 0.2;                                                                   % profile shape factor
x0 = [y0 xc w A mu];                                                        % initial values for peak fitting 
R2tolerance = 0.5;                                                          % minimum value allowed for goodness of fit R2

dmin = 3;                                                                   % minimum number of fitted d values to allow for linear fitting over sin2(psi)


%% Start of data processing
tic                                                                         % tracking processing time

fileinfo = dir(append(dir_data,'*.dat'));                                   % searching for all .dat files from directory
filenames = {fileinfo.name}';                                               % matrix containing all .dat files from directory

psi = [0:10:90]';                                                           % azimuthal angle used for the measurements
sin2psi = sind(psi).^2;                                                     % calculating the corresponding values for sin2(psi)

peak_fit_data = cell(10,1);                                                 % initialise variable
data = cell(10,1);                                                          % initialise variable


%% Reading 10 DPDAK .dat data files into cell for psi = 0 to 90°
for i = 1:10
    data{i} = func_read_dpdak_cake_file(dir_data,filenames{i},nrofscans_total);
end    


%% Fitting all 2theta, intensity (x,y) peak files with Pseudo-Voigt function
lengthx0 = length(x0);
for i = 1:10                                                            	% loop covering all 10 .dat files
    cur_rsquare = zeros(nrofscans,1);
    for j = 1:nrofscans                                                   	% loop covering all scans
        [cur_fitted_curve,cur_gof] = func_pseudo_voigt_peak_fit(data{i}(:,1),data{i}(:,j+1),x0);
                                                                            % fitting data with Pseudo-Voigt function, get fit parameters and goodness of fit
        if cur_gof.rsquare > R2tolerance                                   	% exclude fit results with R2 < R2tolerance
            peak_fit_data{i}(j,:) = coeffvalues(cur_fitted_curve);        	% get fit parameters for current fit                             
            cur_rsquare(j) = cur_gof.rsquare;                              	% get R2 value for current fit = goodness of fit     
        else
            peak_fit_data{i}(j,:) = NaN(1,lengthx0);                       	% replace bad fit parameters with NaN                            
            cur_rsquare(j) = cur_gof.rsquare;                              	% get R2 value for current fit      
        end
    end
    peak_fit_data{i} = [peak_fit_data{i},cur_rsquare];                      % write fit parameters for all scans
end


%% Calculate lattice spacing from fitted peak positions
D = zeros(10,nrofscans);                                                    % matrix for all lattice spacing (d) values
rsquare_peakfit = zeros(10,nrofscans);                                      % matrix for all R2 values

for i = 1:10
    cur_pos = peak_fit_data{i}(:,2);                                        % get fitted peak positions xc
    cur_d = lambda_sync./(2*sind(cur_pos(:,1)./2));                         % calculate lattice spacing from Bragg equation
    peak_fit_data{i} = [peak_fit_data{i},cur_d];                            % add lattice spacing to fit parameter matrix
    
    D(i,:) = peak_fit_data{i}(:,lengthx0+2)';                               % also write lattice spacing in dedicated matrix
    rsquare_peakfit(i,:) = peak_fit_data{i}(:,lengthx0+1)';                 % also write R2 value in dedicated matrix
end


%% Get indexes where peak fit is above R2-tolerance and include in d - sin2(psi) plot for linear fitting%%
NaNidx = ~isnan(D);


%% Linear fit each column of D matrix and return slope, intersection with y-axis, and R2
[k,intersect,rsquare_dfit] = func_linear_fit_sin2psi(sin2psi,D,NaNidx,nrofscans,dmin);           


%% Calculate d0 corresponding to stress free direction using linear fit of d values
d0 = k*sind(psi_stress_free)^2 + intersect;


%% Calculate lattice strain
strain = zeros(10,nrofscans);
for i=1:nrofscans
   strain(:,i) = (D(:,i)-d0(i))./d0(i); 
end    


%% Reshape variables to fit n by m array of scanned points
d0 = reshape(d0,[scan_rows,scan_columns]);                                	% reshape d0 column vector into n by m matrix
k = reshape(k,[scan_rows,scan_columns]);                                   	% reshape slope column vector into n by m matrix
intersect = reshape(intersect,[scan_rows,scan_columns]);                   	% reshape intersect column vector into n by m matrix
rsquare_dfit = reshape(rsquare_dfit,[scan_rows,scan_columns]);             	% reshape rsquare column vector into n by m matrix


%% Calculate in-plane stress for each data point
sigma = k./(half_s2.*d0);                                                  	% calculate stress matrix from slope matrix


%% Export result matrixes
writematrix(sigma,append(dir_data,'sigma.txt'),'Delimiter',';');           	% export of stress matrix to .txt file
writematrix(rsquare_dfit,append(dir_data,'rsquare_dfit.txt'),'Delimiter',';');  
                                                                            % export of R2 matrix from lattice spacing fit to .txt file


%% Plottin the results
% only works for n x m array with n>2 and m>2 in current form

X = 0:scan_stepsize:scan_stepsize*(scan_columns-1);                        	% n horizontal x values spaced by scan stepsize
Y = 0:scan_stepsize:scan_stepsize*(scan_rows-1);                          	% m vertical x values spaced by scan stepsize

f1 = figure;
contourf(X,Y,sigma,1000,'edgecolor','none');                              	% contour plot of stress
set(gca, 'YDir','reverse')
colorbar;                                           
caxis([-4 0])

f2 = figure;
contourf(X,Y,rsquare_dfit,'edgecolor','none');                             	% contour plot of rsquare
colorbar;                                           
set(gca, 'YDir','reverse')
caxis([0 1])

autoArrangeFigures(1);                                                      % function for automated figure window alignment
                                                                            % freely available at Mathworks

toc