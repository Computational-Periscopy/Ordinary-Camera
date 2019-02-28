%% Scene reconstruction using shift and combine method (Fig4 column e)
%  And true occluder positions.

% ------------------------------Pseudo-code------------------------------
% 1. Simulate light transport matrix A given scene geometry parameters.
%    
% 2. Solve TV regularized optimization problem to reconstruct scene
% (FISTA) for 49 postulated occluder shifts, then combine using alpha
% trimmed mean.
% -----------------------------------------------------------------------

% Last Modified by Charles Saunders at Boston University
% 09-Nov-2018 (Clean-up and commented for sharing)

% Manuscript:
%   Saunders, C. and Murray-Bruce, J and Goyal, V.K., 'Computational
%               Periscopy with and Ordinary Digital Camera', Nature, 2018.

% Functions
addpath('Functions')

clear variables;
clc;

TestLetter = 'D11';
numPixels = 1008;

% Parameters
Ndiscr_mon = 4; %Discretization of each scene patch
viewAngleCorrection = 1;
useEstimatedOccPos = false;   %Use estimated occluder position or not
downsamp_factor = 3; %Downsampling of measurements 2^downsamp_factor


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
Mon_xdiscr = linspace(Monitor_xlim(1),Monitor_xlim(2),NumBlocks_col);
Mon_zdiscr = linspace(Monitor_zlim(2),Monitor_zlim(1),NumBlocks_row);

wallparam.wall_matr = wall_matr;
wallparam.wall_point = wall_point;
wallparam.wall_vector_1 = wall_vector_1;
wallparam.wall_vector_2 = wall_vector_2;
wallparam.wall_normal = wall_normal;
wallparam.walln_points = walln_points;

numPixels = floor(numPixels/(2^downsamp_factor));

%% Load data for all test scenes
[test_image1,ground_truth1]=load_image1('image_test_mushroom20.mat',calibParams.filepath,downsamp_factor);   
[test_image2,ground_truth2]=load_image1('image_test_smilehat20.mat',calibParams.filepath,downsamp_factor);
[test_image3,ground_truth3]=load_image1('image_test_colbar20.mat',calibParams.filepath,downsamp_factor);
[test_image4,ground_truth4]=load_image1('image_test_bur20.mat',calibParams.filepath,downsamp_factor);


%%
count = 1;
true_corner = Occ_LLcorner;

scale = D/(D-Occ_LLcorner(2));
px_shiftx = IlluminationBlock_Size(1)/scale;
px_shifty = IlluminationBlock_Size(2)/scale;

%4 pixels in steps of 2
xpos = [-6*px_shiftx, -4*px_shiftx,-2*px_shiftx, 0,  2*px_shiftx, 4*px_shiftx, 6*px_shiftx]; %Offsets of occluder to use
zpos = [-6*px_shifty, -4*px_shifty,-2*px_shifty, 0,  2*px_shifty, 4*px_shifty, 6*px_shifty];


pixel_shift_x = (xpos*scale)/IlluminationBlock_Size(1);
pixel_shift_y = (zpos*scale)/IlluminationBlock_Size(2);

pad = round([max(pixel_shift_x), max(pixel_shift_y)]);

tv_reg_param = 1.3*[0.05,0.05,0.05];

for zoff = 1:length(zpos)
for xoff = 1:length(xpos)
    %%%%%%%%%%%%%

    Occ_LLcorner = true_corner + [xpos(xoff), 0, zpos(zoff)];
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

    tic
    [simA] = simulate_A(wallparam, occ_corner,simuParams, Mon_xdiscr,Mon_zdiscr, 0);
    
    sr = 0.485*0.75*(125*125/((Ndiscr_mon^2)*prod(subblocksperaxis))); 
    sg = 0.485*0.83*(125*125/((Ndiscr_mon^2)*prod(subblocksperaxis))); 
    sb = 0.485*0.89*(125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)));

    final_im1(:,:,:,count) = reconstruct_tv(sr*simA,sg*simA,sb*simA, test_image1,  3.25*tv_reg_param, NumBlocks_sim);
    
    sr =  0.525*0.58*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sg =  0.525*0.72*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sb =  0.525*1.07*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis));
    
    final_im2(:,:,:,count) = reconstruct_tv(sr*simA,sg*simA,sb*simA, test_image2, 3.25*tv_reg_param, NumBlocks_sim);
    
    sr = 0.9*0.75*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sg = 0.9*0.75*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sb = 0.9*0.75*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis));
    
    final_im3(:,:,:,count) = reconstruct_tv(sr*simA,sg*simA,sb*simA, test_image3, 3.5*tv_reg_param, NumBlocks_sim);
    
    sr = 0.65*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sg = 0.85*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis)); 
    sb = 0.95*125*125/((Ndiscr_mon^2)*prod(subblocksperaxis));
    
    final_im4(:,:,:,count) = reconstruct_tv(sr*simA,sg*simA,sb*simA,test_image4, 3.25*tv_reg_param, NumBlocks_sim);
    
    
    final_im1_t(:,:,:,count) = pad_translate(final_im1(:,:,:,count),[pixel_shift_x(xoff),pixel_shift_y(zoff)],pad);
    final_im2_t(:,:,:,count) = pad_translate(final_im2(:,:,:,count),[pixel_shift_x(xoff),pixel_shift_y(zoff)],pad);
    final_im3_t(:,:,:,count) = pad_translate(final_im3(:,:,:,count),[pixel_shift_x(xoff),pixel_shift_y(zoff)],pad);
    final_im4_t(:,:,:,count) = pad_translate(final_im4(:,:,:,count),[pixel_shift_x(xoff),pixel_shift_y(zoff)],pad);

    count = count + 1
    toc
end
end

%%
for i=1:49
    subplot(7,7,i)
    imshow(cat(3,final_im3_t(:,:,1,i),final_im3_t(:,:,2,i),final_im3_t(:,:,3,i)))
end

%%
figure()
subplot(4,4,1)
imshow(ground_truth1)
title('Ground truth')

subplot(4,4,2)
imshow(test_image1./10000000)
title('Measurement')

subplot(4,4,3)
imm1 = cat(3,stack_combine(squeeze(final_im1_t(:,:,1,:))),stack_combine(squeeze(final_im1_t(:,:,2,:))),stack_combine(squeeze(final_im1_t(:,:,3,:))));

crop_coords = [7,6];
crop_size = [29,36];

imshow(imm1);
title('Combined reconstruction')

subplot(4,4,4)
crop_coords = [6,6];
crop_size = [29,36];
imshow(imm1(crop_coords(1):crop_coords(1)+crop_size(1),crop_coords(2):crop_coords(2)+crop_size(2),:));
title('Cropped reconstruction')

%2
subplot(4,4,5)
imshow(ground_truth2)
title('Ground truth')

subplot(4,4,6)
imshow(test_image2./3000000)
title('Measurement')

subplot(4,4,7)
imm2 = cat(3,stack_combine(squeeze(final_im2_t(:,:,1,:))),stack_combine(squeeze(final_im2_t(:,:,2,:))),stack_combine(squeeze(final_im2_t(:,:,3,:))));

imshow(imm2);
title('Combined reconstruction')

subplot(4,4,8)

crop_coords = [6,6];
crop_size = [29,36];

imshow(imm2(crop_coords(1):crop_coords(1)+crop_size(1),crop_coords(2):crop_coords(2)+crop_size(2),:));
title('Cropped reconstruction')

% 3

subplot(4,4,9)
imshow(ground_truth3)
title('Ground truth')

subplot(4,4,10)
imshow(test_image3./3000000)
title('Measurement')

subplot(4,4,11)
imm3 = cat(3,stack_combine(squeeze(final_im3_t(:,:,1,:))),stack_combine(squeeze(final_im3_t(:,:,2,:))),stack_combine(squeeze(final_im3_t(:,:,3,:))));

imshow(imm3)
title('Combined reconstruction')

subplot(4,4,12)
crop_coords = [6,6];
crop_size = [29,36];
imshow(imm3(crop_coords(1):crop_coords(1)+crop_size(1),crop_coords(2):crop_coords(2)+crop_size(2),:)./max(imm3(:)),[]);
title('Cropped reconstruction')

% 4
subplot(4,4,13)
imshow(ground_truth4)
title('Ground truth')

subplot(4,4,14)
imshow(test_image4./3000000)
title('Measurement')

subplot(4,4,15)
imm4 = cat(3,stack_combine(squeeze(final_im4_t(:,:,1,:))),stack_combine(squeeze(final_im4_t(:,:,2,:))),stack_combine(squeeze(final_im4_t(:,:,3,:))));

imshow(imm4);
title('Combined reconstruction')

subplot(4,4,16)
crop_coords = [6,6];
crop_size = [29,36];
imshow(imm4(crop_coords(1):crop_coords(1)+crop_size(1),crop_coords(2):crop_coords(2)+crop_size(2),:));
title('Cropped reconstruction')
