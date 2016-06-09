function [output] = mix(source, i)
%OUTPUT Gives the input source filtered with A(i, j) vector from the A
%matrix given in Ozerov's website

load('sources/mixing_filters_ozerov.mat'); %Ozerov's values
vect1 = squeeze(A(1,i,:));
vect2 = squeeze(A(2,i,:));

output = [filter(vect1, 1, squeeze(source(:, 1))) filter(vect2, 1, squeeze(source(:, 1)))];

end

