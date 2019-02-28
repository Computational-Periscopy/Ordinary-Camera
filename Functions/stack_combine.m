function [final_image] = stack_combine(imagestack)

final_image = zeros(size(imagestack,1),size(imagestack,2));

for ir = 1:size(imagestack,1)
for ic = 1:size(imagestack,2)
   
%     measurements = squeeze(imagestack(ir,ic,:));
%     measurements = measurements(~isnan(measurements));

    
% Distance from median variance
%     med_dist = abs(measurements-median(measurements));
%     selected = measurements(med_dist<mean(med_dist)+std(med_dist)/4);
%     final_image(ir,ic,1) = mean(selected); 

% Percentile
%      measurements = squeeze(imagestack(ir,ic,:));
%      sorted = sort(measurements(~isnan(measurements)));
%      
%      pc = ceil(length(sorted)*0.25); %Percentage to use
%      selected = sorted(floor(length(sorted)/2)-pc:ceil(length(sorted)/2)+pc);
%      
%      final_image(ir,ic,1) = mean(selected);
%      

%  Cluster in range
    med = nanmedian(squeeze(imagestack(ir,ic,:)),1);
    final_image(ir,ic,1) = mean(imagestack(ir,ic,abs(med-imagestack(ir,ic,:))<=0.25));


%     figure()
%     hist(squeeze(imagestack(ir,ic,:)),50)
%     hold on
%     stem(med,1.25,'r')
%     stem(med-0.25,1.25,'r')
%     stem(med+0.25,1.25,'r')
    
end
end
    
end

