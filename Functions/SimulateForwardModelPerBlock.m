function [ Image, MonPattern, Wall_discr ] = SimulateForwardModelPerBlock(NumBlocks,...
    Ndiscr_mon, ActiveBlocks_coord, numPixels, Dist_MonWall, FOV, Occluder,...
    viewAngleCorrection,IlluminationBlock_size, Mon_Offset)
%MONITORSIMULATEFORWARDMODEL simulates the measurements on the visible wall
% (i.e. the camera measurements) for a specified active block (scene patch)
% on the non-visible monitor.
%   Usage:
%   [ Image, MonPattern, Wall_discr ] = MonitorSimulateForwardModel(NumBlocks,
%    Ndiscr_mon, ActiveBlocks_coord, numPixels, Dist_MonWall, FOV, Occluder,
%    viewAngleCorrection,IlluminationBlock_size, Mon_Offset)
%   Input:
%       * NumBlocks (2-vect): Num of possible scene/monitor blocks for
%                   display. (Scene patch discretization essentially).
%       * Ndiscr_mon (scalar): Num of point sources per dimension (of
%                   monitor block). This assumes that each scene block is
%                   Ndiscr_mon*Ndiscr_mon point sources on a uniform grid.
%       * ActiveBlocks_coord (2-vect): [row col] of scene/monitor patch.
%       * numPixels (scalar): Number of camera pixels per color channel to 
%                   simulate. Number of rows of forward matrix (A).
%       * Dist_MonWall (scalar): Distance between monitor and visible wall.
%       * FOV [m] (2-vect): Camera FoV on the visible wall ([x z])
%       * Occluder [m] (3-vector): Position of the lower left corner of the
%                   occluder ([x y z]). 
%       * viewAngleCorrection (boolean): View angle correction for monitor [On/Off].
%       * IlluminationBlock_size[m] (2-vector): Size of monitor blocks ([x,z]).
%       * Mon_Offset[m] (2-vector): Lower left corner of screen ([x, z] only).
%   Output:
%       * Image: Simulated camera measurement (or one column of measurement matrix).
%       * MonPattern: Pattern on the monitor (visualizes ActiveBlocks_coord).
%       * Wall_discr: The [x z] dixretization of visible wall.

% Last Modified by $John Murray-Bruce$ at Boston University.
% v1.0: 20-Nov-2017 9:13:22
% v2.0: 07-Nov-2017 4:46:16 (Clean-up and commented for sharing -JMB)

Occ_present = true;
if isempty(Occluder)
    Occ_present = false;
elseif Occluder==false
    Occ_present = false;
else
    occ_stickpresent = true;
    occstick(:,1) = [-0.0035 0.0035]' + 0.5*(Occluder(1,1)+Occluder(2,1));
    occstick(:,3) = [0 Occluder(1,3)]';
    occstick(:,2) = Occluder(2,2);
end


% DISCRETIZATION
NumBlocks_col = NumBlocks(2);
NumBlocks_row = NumBlocks(1);


% WALL DEFINITION
LambertianWall_xlim = FOV(:,1); %[0.5335 0.988];
LambertianWall_y = Dist_MonWall;
LambertianWall_zlim = FOV(:,2); %[0 0.51];
Wall_xdiscr = linspace(LambertianWall_xlim(1),LambertianWall_xlim(2),numPixels);
Wall_zdiscr = linspace(LambertianWall_zlim(1),LambertianWall_zlim(2),numPixels);


% MONITOR DEFINITION
if nargin<9
    ScreenSize = [0.337 0.271];
    IlluminationBlock_size =0.5*(ScreenSize(1)/(NumBlocks_col+0.5)+...
        (ScreenSize(2)-0.074)/NumBlocks_row);
    IlluminationBlock_size = [IlluminationBlock_size IlluminationBlock_size];
    
    MonitorOffset_x = 0.091;
    MonitorOffset_z = 0.1445;
elseif nargin==9
    MonitorOffset_x = 0.091;
    MonitorOffset_z = 0.1445;
elseif nargin==10
    MonitorOffset_x = Mon_Offset(1);
    MonitorOffset_z = Mon_Offset(2);
end


Monitor_xlim = [0 NumBlocks_col]*IlluminationBlock_size(1) + MonitorOffset_x;
Monitor_y = 0;
Monitor_zlim = [0 NumBlocks_row]*IlluminationBlock_size(2) + MonitorOffset_z;
N_monDiscr_row = Ndiscr_mon*NumBlocks_row;
N_monDiscr_col = Ndiscr_mon*NumBlocks_col;
Mon_xdiscr = linspace(Monitor_xlim(1),Monitor_xlim(2),N_monDiscr_col);
Mon_zdiscr = linspace(Monitor_zlim(2),Monitor_zlim(1),N_monDiscr_row);


% INITIALIZE ALL DESIRED SCENE/MONITOR PATCHES TO BE COMPUTED (number of columns of A matrix)
ActiveBlocks_matrix = zeros(NumBlocks_row,NumBlocks_col);
ActiveBlocks_matrix(ActiveBlocks_coord(:,1),ActiveBlocks_coord(:,2)) = 1; % Intensity of blocks.
MonitorPattern_comp = kron(ActiveBlocks_matrix,ones(Ndiscr_mon,Ndiscr_mon));
MonPattern = MonitorPattern_comp;
MonitorPattern_comp = (fliplr(MonitorPattern_comp))';
MonitorPattern_comp = MonitorPattern_comp(:);

%Discretize visible wall (within the camera FOV)
[XX_wall,YY_wall,ZZ_wall] = meshgrid(Wall_xdiscr,LambertianWall_y,Wall_zdiscr);
Pos_wall = [XX_wall(:),YY_wall(:),ZZ_wall(:)];
%Discretize hidden-scene (with the LCD screen visible region)
[XX_mon,YY_mon,ZZ_mon] = meshgrid(Mon_xdiscr,Monitor_y,Mon_zdiscr);
Pos_mon = [XX_mon(:),YY_mon(:),ZZ_mon(:)];

Normal_monitor = [0, 1, 0];
Normal_wall = [0, -1, 0];

% Preallocations
tempOnes = ones(numPixels*numPixels,1);
LL = zeros(numPixels*numPixels,1);

% Compute the A matrix
for ii=1:N_monDiscr_col*N_monDiscr_row
    if MonitorPattern_comp(ii)>0
        tDiff =(tempOnes*Pos_mon(ii,:) - Pos_wall);
        temp1 = (tDiff*Normal_wall')./(sum(tDiff.^2,2));
        temp1(temp1<0)= 0;
        temp2 = (-tDiff*Normal_monitor')./(sum(tDiff.^2,2).^2);
        temp2(temp2<0)= 0;
        
        % If occluder is present incorporate the visibility of occluder
        if Occ_present==true
            % Visibility along x-dimension
            [Vx] = ComputeOccluderShadow(Pos_mon(ii,1),Wall_xdiscr,Dist_MonWall,Occluder(1,2),Occluder(:,1));
            
            % Visibility along z-dimension
            [Vz] = ComputeOccluderShadow(Pos_mon(ii,3),Wall_zdiscr,Dist_MonWall,Occluder(2,2),Occluder(:,3));
            
            % Outter product below yields 2D visibility matrix (fast)
            visibilityMat = (1-(1-Vz)'*(1-Vx))';
            visibilityMat = visibilityMat(:);
            
            % Occluder is suspended using a stand, so compute the stand's
            % vibility also and combine with above rectangular visibility.
            if occ_stickpresent==true
                [Vx_s] = ComputeOccluderShadow(Pos_mon(ii,1),Wall_xdiscr,Dist_MonWall,occstick(1,2),occstick(:,1));
                [Vz_s] = ComputeOccluderShadow(Pos_mon(ii,3),Wall_zdiscr,Dist_MonWall,occstick(2,2),occstick(:,3));
                visibilityMat = (1- (((1-Vz)'*(1-Vx))|((1-Vz_s)'*(1-Vx_s))))';
                visibilityMat = visibilityMat(:);
            end
        else
            visibilityMat=1;
        end
        
        % View angle correction to model LCD view angle variations
        if viewAngleCorrection==true
            MM = ViewingAngleFactor(Pos_mon(ii,:),FOV,Dist_MonWall,numPixels);
        else
            MM=1;
        end
        LL(:,1) = LL(:,1) + temp1.*temp2.*visibilityMat*MonitorPattern_comp(ii).*(MM(:));
    end
end

Image = reshape(sum(LL,2),numPixels,numPixels).';
Wall_discr = [Wall_xdiscr;Wall_zdiscr];

end


function [M] = ViewingAngleFactor(MonitorPixel_xyz, Camera_FOV, D, numPixels)
%Models the view angle variations for the LCD monitor.
powexp = 20;
FOV_zdiscr = linspace(Camera_FOV(2,2),Camera_FOV(1,2),numPixels)';
FOV_xdiscr = linspace(Camera_FOV(2,1),Camera_FOV(1,1),numPixels)';
Mz = cos(1*atan((MonitorPixel_xyz(3)-FOV_zdiscr)./(D))).^powexp;
Mx = cos(atan((MonitorPixel_xyz(1)-FOV_xdiscr)./(D))).^0;
M = (Mx*Mz');
end


function [Vx,X1,X2] = ComputeOccluderShadow(x_locs,xx,D_wall,Occl_d,Occ_edges)
%ComputeOccluderShadow computes the visibility function of a planar
%occluder defined by the occluder distance from monitor 'Occl_d', and edges
% i.e. Occ_edges, from the point of view of x_locs. This is a 1D
% computation. An outter product of two 1D visibilities gives a full 2D
% visibility matrix. Separability is exploited for speed.
%
%   Usage:
%       [Vx,X1,X2] = ComputeOccluderVisibility(x_locs,xx,D_wall,Occl_d,Occ_edges)
%
%   Input:
%       * x_locs:		Position of hidden-scene patch in 2D.
%       * xx:           Discretization of visible wall (i.e. Camera FOV).
%       * D_wall:       Distance between visible wall and monitor (or scene plane).
%       * Occ_d:        Distance.
%       * Occ_edges:    Edges of the occluder.
%   Output:
%       * Vx:           Visibility matrix of the occluder on the visible wall.

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 28-Aug-2017 11:13:11
% v2.0 07-Nov-2018 17:04:43

xx = xx(:).';              % Enforce that xx is a row vector.
N = length(xx);
[Ns] = length(x_locs);
XX = repmat(xx,Ns,1);

X1 = (D_wall*(Occ_edges(1) - x_locs)./Occl_d) +  x_locs;
X1 = repmat(X1,1,N);
X2 = (D_wall*(Occ_edges(2) - x_locs)./Occl_d) +  x_locs;
X2 = repmat(X2,1,N);

Vx = ones(Ns,N);
Vx( (XX>=X1)&(XX<=X2) ) =0;



end

