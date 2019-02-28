function [output] = pad_translate(image,shift,pad)

impad = padarray(image,[pad(1),pad(2),0],0,'both');
output = imtranslate(impad,shift);

end

