%% This script details measurements configurations for the experimental tests
% performed.
%

filename = 'TestPosD11';
NumBlocks_sim = [29 36];    % Number of scene patches [hor, vert]
D = 1.03;

% OCCLUDER DEFINITION
Occ_size = [0.077 0 0.075]; % Size of the flat occluder (x,y,z)

% Please note here that the oocluder's y-position is stated as (D-p_y), where
% D = 1.03 (D: distance between hidden-scene and visible wall
% planes) and p_y is the distance between the visible wall plane and
% occluder plane.
if useEstimatedOccPos
    % Estimated occluder positions for the different scenes
    %             Occ_LLcorner = [0.4583    D-0.4892    0.2026];  % Mushroom
    %             Occ_LLcorner = [0.4569    D-0.4556    0.2080];  % Tommy
    %             Occ_LLcorner = [0.4733    D-0.4639    0.2072];  % BU red
    %             Occ_LLcorner = [0.4693    D-0.4674    0.2080];  % RGB
else
    Occ_LLcorner = [0.470 D-0.460 0.2040];  % True/measured occluder position
end
% Define bottom-right and top-left edges of occluder
Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

%CAMERA CONFIG
FOV_size = [0.4372 0.4372];
FOV_LLCorner = [0.5128 0.0482];
FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

% MONITOR CONFIG
Mon_Offset = [0.0210 0.136];        % (x,z)-position of the lower-right edge of usuable portion of LCD screen

NumBlocks_col = NumBlocks_sim(2);   % Number of scene patches (horizontally)
NumBlocks_row = NumBlocks_sim(1);   % Number of scene patches (vertically)
ScreenSize = [0.408 0.3085];        % Size of the LCD display [m]
ScreenResolution = [1280 1024];     % [hor, ver] pixels
NumPixelsPerMonitorBlock = 35;      % Number of monitor pixels in each hidden-scene patch
PixelSize_m = (ScreenSize./ScreenResolution);   % Dimension of each scene patch [m]

Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_sim(1));
Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_sim(2));

IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock;



% Set parameters for simulating the forward model (i.e. A matrix)
simuParams.NumBlocks = NumBlocks_sim;
simuParams.Ndiscr_mon = Ndiscr_mon;
simuParams.numPixels = numPixels;
simuParams.D = D;
simuParams.FOV_cord = FOV_cord;
simuParams.Occluder = Occluder;
simuParams.viewAngleCorrection = viewAngleCorrection;
simuParams.IlluminationBlock_Size = IlluminationBlock_Size;
simuParams.Mon_Offset = Mon_Offset;


if ismac
    datafilepath = './Data/TestPosD11/';
elseif ispc
    datafilepath = '.\Data\TestPosD11\';
end
