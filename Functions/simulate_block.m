%% Simulates one monitor block. GPU accelerated.
%% Charles Saunders at Boston University

function [image] = simulate_block(wall_mat,wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos, occ_corner,MM)


    % Squares each element in input
    function [vecsqr] = calc(vec)
        vecsqr =  vec.^2;
    end

    % Vector from light position to wall position.
    vec = repelem(wall_mat([1,3],:),1,1,size(light_source_pos,1)) - repelem(reshape(light_source_pos(:,[1,3])',[2,1,length(light_source_pos)]),1,walln_points,1);
    % Distance squared
    vs = arrayfun(@calc, vec);
    
    % y component is constant, so:
    vs2 = ones(1,walln_points)*((wall_mat(2,1)-light_source_pos(1,2)).^2);

    % Calculate distance squared again (to normalize dot products and for
    % distance losses).
    tmp = arrayfun(@calc,repmat(vs(1,:,:) + vs2,[walln_points,1,1]) + repmat(permute(vs(2,:,:),[2,1,3]),[1,walln_points,1]));
    
    % Final intensity
    intensity = MM.*repmat(vs2,[walln_points,1,length(light_source_pos)])./tmp;
    
    % Calculate occluder positions
    image = zeros(walln_points,walln_points,'gpuArray');
    if size(occ_corner)>0
    for i = 1:size(vec,3)
        occluder_image = simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos(i,:), occ_corner(:,:,1));
        for o = 2:size(occ_corner,3)
            occluder_image = occluder_image.*simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos(i,:), occ_corner(:,:,o));
        end
        
        image = image + intensity(:,:,i).*occluder_image;
    end
    else
        image = sum(intensity,3);
    end

end

