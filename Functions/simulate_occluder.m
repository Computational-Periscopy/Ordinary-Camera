%% Generates occluder mask
%% Charles Saunders at Boston University

function [ occluder_image ] = simulate_occluder(wall_point, wall_vector_1, wall_vector_2, wall_normal, walln_points, light_source_pos, occ_corner)

    % Find projection of occluder corners onto wall
    occ_corner_proj = zeros(size(occ_corner,1),3);
    for i = 1:size(occ_corner,1)
        occ_corner_proj(i,:) = ray_plane_intersect_y(light_source_pos,occ_corner(i,:),wall_point,wall_normal);
    end

    % Coordinates in terms of the FOV
    occ_image_coords_1 = (occ_corner_proj(:,1) - (wall_point(1) - wall_vector_1(1))) * walln_points/(wall_vector_1(1)*2);
    occ_image_coords_2 = (occ_corner_proj(:,3) - (wall_point(3) - wall_vector_2(3))) * walln_points/(wall_vector_2(3)*2);
    
    % Generate mask
     k = convhull(occ_image_coords_1,occ_image_coords_2,'simplify', true);
     occluder_image = ~poly2mask(occ_image_coords_1(k),occ_image_coords_2(k),walln_points,walln_points);

    
end