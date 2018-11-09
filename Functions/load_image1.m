function [image,ground_truth] = load_image1(file,path,downsampling)
% Loads an image, removes background and splits into colour channels


iml = load([path,file]);
imload = double(iml.image);
ground_truth = iml.ground_truth;
measure_shut = get_shutter_speed(iml);



if ~exist('calshut','var')
    calshut = measure_shut;
end

imcol = get_color_image(imload);

testimr = (calshut/measure_shut)*downsample_2(imcol(:,:,1),downsampling);
testimg = (calshut/measure_shut)*downsample_2(imcol(:,:,2),downsampling);
testimb = (calshut/measure_shut)*downsample_2(imcol(:,:,3),downsampling);

image(:,:,1) = testimr;
image(:,:,2) = testimg;
image(:,:,3) = testimb;
end

function [sp] = get_shutter_speed(iml)

speed = iml.shutter_speed;

if isfield(iml,'averaged_iterations')
    sp = speed*iml.averaged_iterations;
else
    sp = speed*iml.iterations;
    
end
end

function [col_image] = get_color_image(im)
	red = im(1:2:end,1:2:end);
    green = (im(2:2:end,1:2:end)+im(1:2:end,2:2:end))/2;
    blue = im(2:2:end,2:2:end);
    
    col_image(:,:,1) = red;
    col_image(:,:,2) = green;
    col_image(:,:,3) = blue;
end

function [image] = downsample_2(im,iterations)
    %Downsampling 2^iterations times without averaging

    for i=1:iterations
        im1 = im(1:2:end,1:2:end);
        im2 = im(1:2:end,2:2:end);
        im3 = im(2:2:end,1:2:end);
        im4 = im(2:2:end,2:2:end);

        image = (im1+im2+im3+im4)/4;
        im = image;
    end
    
    if iterations == 0 
        image = im;
    end
end
