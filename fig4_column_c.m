%% Scene reconstruction using differencing method (Fig4 column c)
%
% ------------------------------Pseudo-code------------------------------
% 1. Simulate light transport matrix A given scene geometry parameters.
%    
% 2. Solve TV regularized optimization problem to reconstruct scene
% (FISTA).
% -----------------------------------------------------------------------

% Last Modified by Charles Saunders at Boston University
% 09-Nov-2018 (Clean-up and commented for sharing)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.

%%
% Functions
addpath('Functions')

clear variables;

TestLetter = 'D11'; % Picks the test data to use
numPixels = 1008; % Number of pixels in camera measurement

% Parameters
Ndiscr_mon = 6; %Discretization of each scene patch
downsamp_factor = 3; %Downsampling of measurements 2^downsamp_factor
viewAngleCorrection = 1;  %True/False
useEstimatedOccPos = true; %Use estimated occluder position or not

load_experiment_config_data_reconstruction

% Data path
if ismac
    calibParams.filepath = './Data/TestPosD11/';
elseif ispc
    calibParams.filepath = '.\Data\TestPosD11\';
end
%%%%%%% Setup %%%%%%%% (Nothing to manually set here)

% Wall/imaging plane
wall_point = [FOV_LLCorner(1) + FOV_size(1)/2,D,FOV_LLCorner(2)+FOV_size(2)/2];  %Point on plane
wall_vector_1 = [FOV_size(1)/2,0,0]; %Vector defining one direction of FOV (and extent)
wall_vector_2 = [0,0,FOV_size(2)/2]; %Vector defining the orthogonal direction (and extent)
wall_normal = cross(wall_vector_1,wall_vector_2);
wall_normal = wall_normal./norm(wall_normal);

walln_points = floor(numPixels/(2^downsamp_factor)); %Number of points to render in each direction

% Discretize imaging plane
f_imageplane = gpuArray((zeros(walln_points)));
wall_vec = (-1:2/(walln_points-1):1);
wall_matr(1,:) = gpuArray((wall_point(1) + wall_vec*wall_vector_1(1) + wall_vec*wall_vector_2(1)));
wall_matr(2,:) = gpuArray((wall_point(2) + wall_vec*wall_vector_1(2) + wall_vec*wall_vector_2(2)));
wall_matr(3,:) = gpuArray((wall_point(3) + wall_vec*wall_vector_1(3) + wall_vec*wall_vector_2(3)));

Monitor_xlim = [0 NumBlocks_col]*IlluminationBlock_Size(1) + Mon_Offset(1);
Monitor_y = 0;
Monitor_zlim = [0 NumBlocks_row]*IlluminationBlock_Size(2) + Mon_Offset(2);
Mon_xdiscr = (linspace(Monitor_xlim(1),Monitor_xlim(2),NumBlocks_col));
Mon_zdiscr = (linspace(Monitor_zlim(2),Monitor_zlim(1),NumBlocks_row));

wallparam.wall_matr = wall_matr;
wallparam.wall_point = wall_point;
wallparam.wall_vector_1 = wall_vector_1;
wallparam.wall_vector_2 = wall_vector_2;
wallparam.wall_normal = wall_normal;
wallparam.walln_points = walln_points;

numPixels = floor(numPixels/(2^downsamp_factor));

%% Load data

% Available hidden test scenes data choose one [scene name: 'key']
%           * RGB bars scene:       'rgb'
%           * Text 'BU' scene:      'bu'
%           * Mushroom scene:       'mushroom'
%           * Tommy scene:          'tommy'

scene = 'bu';

switch scene
    case 'mushroom'
        [test_image1,ground_truth1]=load_image1('image_test_mushroom20.mat',calibParams.filepath,downsamp_factor);   
        Occ_LLcorner = [0.4583 0.5408 0.2026]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e08 * [0.51    0.561    2.448]; %TV regularization parameter
        
        %Discretization dependent scaling of A matrix
        sr = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.5313e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 16250/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'tommy'
        [test_image1,ground_truth1]=load_image1('image_test_smilehat20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4569 0.5744 0.2080]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e07 * [5.72    6.76    5.72];  %TV regularization parameter
        
         %Discretization dependent scaling of A matrix
        sr = 1.1406e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.3594e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.9063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'bu'
        [test_image1,ground_truth1]=load_image1('image_test_bur20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4733   0.5661    0.2072]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e07*[5   25   25];  %TV regularization parameter
        
         %Discretization dependent scaling of A matrix
        sr = 12500/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 15625/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.7188e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
    case 'rgb'
        [test_image1,ground_truth1]=load_image1('image_test_colbar20.mat',calibParams.filepath,downsamp_factor);
        Occ_LLcorner = [0.4693 0.5629 0.2080]; %Estimated occluder position from localization script
        
        tv_reg_param = 1e06 * [52.5   50   47.5];  %TV regularization parameter
        
        %Discretization dependent scaling of A matrix
        sr = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.9; 
        sg = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*0.98; 
        sb = 1.4063e+04/(Ndiscr_mon^2)*prod(subblocksperaxis)*1.04;
end

%%
% Occluder 
occ_corner(1,:,1) = Occ_LLcorner;
occ_corner(2,:,1) = Occ_LLcorner + [Occ_size(1), 0, 0];
occ_corner(3,:,1) = Occ_LLcorner + [Occ_size(1), 0, Occ_size(3)];
occ_corner(4,:,1) = Occ_LLcorner + [0, 0, Occ_size(3)];

occ_corner(1,:,2) = Occ_LLcorner + [Occ_size(1)/2-0.00275, 0, 0];
occ_corner(2,:,2) = Occ_LLcorner + [Occ_size(1)/2+0.00275, 0, 0];
occ_corner(3,:,2) = Occ_LLcorner + [Occ_size(1)/2-0.00275, 0, -Occ_LLcorner(3)];
occ_corner(4,:,2) = Occ_LLcorner + [Occ_size(1)/2+0.00275, 0, -Occ_LLcorner(3)];
%%%%%%%%%%%%%%


%% Simulate Transport Matrix
disp('Simulating transport matrix...')
[simA] = simulate_A(wallparam, (occ_corner),simuParams, Mon_xdiscr,Mon_zdiscr, 0);

%% Reconstruction
final_im1 = reconstruct_tv_it_cbg(simA,[sr,sg,sb],  test_image1, tv_reg_param, NumBlocks_sim, [0,0,0]);

%% Plots
figure()

subplot(1,2,1)
imshow(ground_truth1/255)
title('Ground truths')
subplot(1,2,2)
imshow(final_im1(:,:,:))
