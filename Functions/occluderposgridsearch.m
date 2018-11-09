function [p_est, p_est_rgb, K_vals] = occluderposgridsearch(meas,simuParams,Occ_size,...
    gridvals,II,ds_factor,sigma_th)
%OCCLUDERPOSITIONGRIDSEARCH finds the position of the hidden occluder by
% performing a search on the positions specified by gridvals. This is
% achieved by maximizing the functional (see Eqn. (2) in the manuscript):
% 
%  (2)   || A(p_o) (A(p_o)^{T} A(p_o))^{-1} A(p_o)^{T} * y ||_{2}^{2}
%
% with respect to p_o (where p_o is each (x,y,z)-position in gridvals.
% See the method section of online manuscript and supplementary material
% for full details of the localization algorithm.
%
%
%   Usage:
%       [p_est, objectivefunc, p_est_rgb] = occluderposgridsearch(meas,simuParams,Occ_size,
%                                               gridvals,II,ds_factor,sigma_th)
%   Input:
%       * meas (struct obj):       Measured data (i.e. the single snapshot taken by camera).
%                                   -meas.r (matrix): red channel data
%                                   -meas.g (matrix): occluder position in 3D-space
%                                   -meas.b (matrix): occluder position in 3D-space
%       * simuParams (struct obj): Containing the parameters used for simulating the forward model A.
%                                   -simuParams.Occluder: occluder position in 3D-space
%                                   -simuParams.NumBlocks: occluder position in 3D-space
%                                   -simuParams.Ndiscr_mon: occluder position in 3D-space
%                                   -simuParams.numPixels: occluder position in 3D-space
%                                   -simuParams.D: occluder position in 3D-space
%                                   -simuParams.viewAngleCorrection: occluder position in 3D-space
%                                   -simuParams.IlluminationBlock_Size: occluder position in 3D-space
%                                   -simuParams.Mon_Offset: LCD position in (x,z)-plane
%       * Occ_size (1 x 2 matrix):  Dimensions of occluder [Length (x-axis), Height (z-axiz)] m.
%       * gridvals (3 x N):         Gridpoints along x, y and z coordinates for candidate occluder positions p_o.
%       * II (1 x 3):               Number of gridpoints along each x, y, z.
%       * ds_factor (scalar):       Amount of downsampling applied onto the camera measurement (=number of rows of A).
%       * sigma_th (scalar):        Normalized singular values threshold used.
%
%   Output:
%       * p_est (1x3 matrix):                   Estimated occluder position.
%       * p_est_rgb (3x3 matrix):               Est. occluder positions per colour channel.
%       * K_vals (II_1 x II_2 x II_3 matrix):	The number of normalized singular vectors
%                                             larger than sigma_th per candidate position.

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 05-Sep-2017 10:28:32
% v2.0 11-Nov-2018 14:09:47 (Clean-up and commented for sharing -JMB)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.






% Gather grid points of (p_o)_x, (p_o)_y and (p_o)_z for the search.
xhatvals = gridvals(1,:);
yhatvals = gridvals(2,:);
zhatvals = gridvals(3,:);

% Preallocate arrays for storing norm (2) for each candidate position p_o.
proj_yOnRangeA_r = zeros(II(1),II(2),II(3));
proj_yOnRangeA_g = zeros(II(1),II(2),II(3));
proj_yOnRangeA_b = zeros(II(1),II(2),II(3));

% Preallocate array for storing number of singular values retained.
K_vals = zeros(II(1),II(2),II(3));

for iix=1:II(1)
    for iiy=1:II(2)
        for iiz=1:II(3)
            % Define candidate occuder position (p_o)
            Occ_LLcorner = [xhatvals(iix) yhatvals(iiy) zhatvals(iiz)];
            Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
            simuParams.Occluder = Occluder;
            
            % Simulate corresponding forward matrix, i.e. A(p_o)
            [ simA, ~ ] = SimulateA_OccluderEstimation(simuParams, ds_factor);
            
            % Compute economical SVD of forward matrix A(p_o)
            [U, S, ~] = svd(simA, 'econ');
            vecS = diag(S)/S(1);            % Normalize singular values by largest
            K = min(sum(vecS>sigma_th),10);	% Number of largest singular values retained
            K_vals(iix,iiy,iiz) = sum(vecS>sigma_th); % Save K values for viewing
            
            % Compute norm of low dimensional projection of measurements,
            % y, along range space of A for each colour channel.
            proj_yOnRangeA_r(iix,iiy,iiz) = norm(U(:,1:K)'*(meas.r(:)));
            proj_yOnRangeA_g(iix,iiy,iiz) = norm(U(:,1:K)'*(meas.g(:)));
            proj_yOnRangeA_b(iix,iiy,iiz) = norm(U(:,1:K)'*(meas.b(:)));
        end
    end
end

% Retrieve the (p_o) value that yields the maximum energy for red channel
[~,Indx] = max(proj_yOnRangeA_r(:));
[I_x, I_y, I_z] = ind2sub(size(proj_yOnRangeA_r),Indx);
p_est_rgb(1,:) = [xhatvals(I_x) yhatvals(I_y) zhatvals(I_z)];

% Retrieve the (p_o) value that yields the maximum energy for green channel
[~,Indx] = max(proj_yOnRangeA_g(:));
[I_x, I_y, I_z] = ind2sub(size(proj_yOnRangeA_g),Indx);
p_est_rgb(2,:) = [xhatvals(I_x) yhatvals(I_y) zhatvals(I_z)];

% Retrieve the (p_o) value that yields the maximum energy for blue channel
[~,Indx] = max(proj_yOnRangeA_b(:));
[I_x, I_y, I_z] = ind2sub(size(proj_yOnRangeA_b),Indx);
p_est_rgb(3,:) = [xhatvals(I_x) yhatvals(I_y) zhatvals(I_z)];

% Take their average
p_est = mean(p_est_rgb);
