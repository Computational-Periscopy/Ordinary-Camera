function [point] = ray_plane_intersect_y(ray_start, ray_end, plane_point, plane_normal)
%% Ray intersection with plane. Optimized for y-axis-alligned planes
%   ray_start = 3d vector, x,y,z
%   ray_end = 3d vector, x,y,z
%   square_coordinates = x,y,z; x,y,z; x,y,z; x,y,z;

% Charles Saunders at Boston University

ray_direction = ray_end-ray_start;
ray_direction = ray_direction./norm(ray_direction,2);

dist = (ray_start(2) - plane_point(2)) / -ray_direction(2);
if dist>0
    point = ray_start + dist * ray_direction;
else  % We don't intersect the plane, but we'll pretend we do and project to get shadows
    point = ray_start - dist * ray_direction;
    point = point + point(2).*plane_normal + plane_point;
end

end

