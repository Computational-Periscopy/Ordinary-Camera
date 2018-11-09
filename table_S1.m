%% Occluder position estimation: Grid search method on cost function
% Starts with a coarse grid and then adaptively refines grid based on most
% recent estimates.
%
% ------------------------------Pseudo-code------------------------------
% 1. Search (px,py,pz) jointly on a coarse 5x5x5 grid. Using the threshold
%    sigma_th1 =0.75, to select the largest K normalized singular values.
% 2. For each coordinate (say px) refine estimate on finer grid and keep
%    others fixed (py, pz). Repeat for other coordinates (in the order px, 
%    pz, py). At the (ii+1)-th iteration, the grid for each coordinate has
%    limits [current_p_est-range_vals(ii), current_p_est+range_vals(ii)].
%    In this way the new search region is centred on the most recent
%    estimate of the occluder position.
% 3. Repeat step 2 with finer grid.
% -----------------------------------------------------------------------

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 23-Jan-2018 10:31:32
% v2.0 08-Nov-2018 14:14:57 (Clean-up and commented for sharing -JMB)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.

% Functions
addpath('Functions')

clear variables;
close all; clc;

% Available hidden test scenes data choose one from [scene name: 'key']
%           * RGB bars scene:       'rgb'
%           * Text 'BU' scene:      'bu'
%           * Mushroom scene:       'mushroom'
%           * Tommy scene:          'tommy'

testscene = 'mushroom';
numPixels = 1008;

% MONITOR DISCRETIZATION
Ndiscr_mon = 1;             % First coarse grid!
viewAngleCorrection = true;
useEstimatedOccPos = false;
load_experiment_config_data_localization;

downsamp_factor = 4;
bgSub = 0;

switch lower(testscene)
    case 'rgb'
        [camera_capture,ground_truth1]=load_image1('image_test_colbar20.mat',datafilepath,downsamp_factor);
        n_iter = 4;
        sigma_th1 = 0.75;
        sigma_th2 = 0.2;
    case 'bu'
        [camera_capture,ground_truth1]=load_image1('image_test_bur20.mat',datafilepath,downsamp_factor);
        n_iter = 3;
        sigma_th1 = 0.75;
        sigma_th2 = 0.5;
    case 'mushroom'
        [camera_capture,ground_truth1]=load_image1('image_test_mushroom20.mat',datafilepath,downsamp_factor);
        n_iter = 3;
        sigma_th1 = 0.75;
        sigma_th2 = 0.2;
    case 'tommy'
        [camera_capture,ground_truth1]=load_image1('image_test_tommy20.mat',datafilepath,downsamp_factor);
        n_iter = 4;
        sigma_th1 = 0.75;
        sigma_th2 = 0.5;
    otherwise
        disp('No such test scene exists, try ''RGB,')
end




meas.r = camera_capture(:,:,1);
meas.g = camera_capture(:,:,2);
meas.b = camera_capture(:,:,3);

N_oneparam = 30;


%
p_est_all = zeros(n_iter,3);
range_vals = [0.025 0.0125 0.00625 0.003125]*4;
for ii=1:n_iter
    a = tic;
    
    if ii<2
        % Global search step
        IIx = 5; IIy = 5; IIz = 5;
        simuParams.Ndiscr_mon = 2;
        xhatvals = linspace(0.1,0.6,IIx);
        yhatvals = linspace(0.1,0.6,IIy);
        zhatvals = linspace(0.0,0.3,IIz);
        II = [IIx IIy IIz];
        gridvals(1,:) =  xhatvals;
        gridvals(2,:) =  yhatvals;
        gridvals(3,:) =  zhatvals;
        [p_est, ~] = occluderposgridsearch(meas,simuParams,Occ_size,...
            gridvals,II,downsamp_factor,sigma_th1);
    else
        % Corordinate descent step
        simuParams.Ndiscr_mon = 1;
        clear gridvals;
        arange_val = range_vals(ii-1);
        sigma_th = sigma_th2/(10^(ii-2));
        
        gridvals(1,:) =  linspace(p_est(1) - arange_val, p_est(1) + arange_val, N_oneparam);
        xhatvals = gridvals(1,:);
        gridvals(2,:) =  p_est(2); gridvals(3,:) =  p_est(3);
        II = [N_oneparam 1 1];
        [p_est, ~] = occluderposgridsearch(meas,simuParams,Occ_size,...
            gridvals,II,downsamp_factor,sigma_th);
        
        gridvals(3,:) =  linspace(p_est(3) - arange_val, p_est(3) + arange_val, N_oneparam);
        zhatvals = gridvals(3,:);
        gridvals(1,:) =  p_est(1); gridvals(2,:) =  p_est(2);
        II = [1 1 N_oneparam];
        [p_est, ~] = occluderposgridsearch(meas,simuParams,Occ_size,...
            gridvals,II,downsamp_factor,sigma_th);
        
        simuParams.Ndiscr_mon = 10;
        gridvals(2,:) =  linspace(p_est(2) - arange_val, p_est(2) + arange_val, N_oneparam);
        yhatvals = gridvals(2,:);
        gridvals(1,:) =  p_est(1); gridvals(3,:) =  p_est(3);
        II = [1 N_oneparam 1];
        [p_est, ~] = occluderposgridsearch(meas,simuParams,Occ_size,...
            gridvals,II,downsamp_factor,sigma_th);
        
        gridvals(1,:) =  xhatvals;
        gridvals(2,:) =  yhatvals;
        gridvals(3,:) =  zhatvals;
    end
    
    p_est_all(ii,:)  = p_est;
    disp(['Interation: ',num2str(ii),', Time Elapsed (s): ' num2str(toc(a))]);
end

% Correct reference point for y-axis (to match manuscript's Fig. 1)
p_est(2) = D - p_est(2);
disp('Estimated Occluder Position [(hat{p}_o)_x, (hat{p}_o)_y, (hat{p}_o)_z] = p_est')

p_est
