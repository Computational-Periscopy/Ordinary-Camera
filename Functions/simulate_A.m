%% Simulate the A matrix using forward model
%% Charles Saunders and John Murray-Bruce at Boston University

function [A] = simulate_A(wallparam, occ_corner, simuParams, Mon_xdiscr, Mon_zdiscr, Monitor_depth)
%% Simulates the full A matrix. GPU accelerated

wall_matr = wallparam.wall_matr;
wall_point = wallparam.wall_point;
wall_vector_1 = wallparam.wall_vector_1;
wall_vector_2 = wallparam.wall_vector_2;
wall_normal = wallparam.wall_normal;
walln_points = wallparam.walln_points;

NumBlocks_row = simuParams.NumBlocks(1);
NumBlocks_col = simuParams.NumBlocks(2);
Ndiscr_mon = simuParams.Ndiscr_mon;

blockcount = 1;

A = (zeros(walln_points^2, NumBlocks_row*NumBlocks_col,'gpuArray'));

for mc = NumBlocks_col:-1:1 %For each scene patch
    for mr = 1:NumBlocks_row
        lightposy = Monitor_depth;

        % View angle model
        if simuParams.viewAngleCorrection==1
            MM = ViewingAngleFactor([Mon_xdiscr(mc)-simuParams.IlluminationBlock_Size(1)/2, lightposy ,Mon_zdiscr(mr)-simuParams.IlluminationBlock_Size(2)/2],wall_matr(1,:), wall_matr(3,end:-1:1), simuParams.D);
        elseif simuParams.viewAngleCorrection==0
            MM = ones(length(wall_matr(1,:)),length(wall_matr(3,end:-1:1)));
        end

        % Discretize a monitor block
        [lightposx, lightposz] = meshgrid( Mon_xdiscr(mc):-simuParams.IlluminationBlock_Size(1)/Ndiscr_mon:Mon_xdiscr(mc)-simuParams.IlluminationBlock_Size(1)+(simuParams.IlluminationBlock_Size(1)/Ndiscr_mon)/10, Mon_zdiscr(mr):-simuParams.IlluminationBlock_Size(2)/Ndiscr_mon:Mon_zdiscr(mr)-simuParams.IlluminationBlock_Size(2)+(simuParams.IlluminationBlock_Size(2)/Ndiscr_mon)/10); 
        
        % Simulate the block
        image = flipud(simulate_block(wall_matr,wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, [lightposx(:),ones(Ndiscr_mon^2,1).*lightposy,lightposz(:)], occ_corner, repmat(MM,[1,1,Ndiscr_mon^2])));

        % Add to A matrix as column
        A(:,blockcount) = image(:);
        blockcount = blockcount + 1;
        
    end
end

A = gather(A);

end

        


function [M] = ViewingAngleFactor(MonitorPixel_xyz, FOV_xdiscr,FOV_zdiscr, D)
powexp = 18;%18;%20; %5;
ang = atan((MonitorPixel_xyz(3)-FOV_zdiscr)./(D));
Mz = cos(ang).^powexp;

Mx = cos(atan((MonitorPixel_xyz(1)-FOV_xdiscr)./(D))).^1;

M = (Mx'*Mz)';
end

