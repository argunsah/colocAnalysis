function [x,y] = convert_1D_to_2D_locs(vals,xsize)

for i = 1:length(vals)
    if vals(i) < xsize
        x(i) = 1;
        y(i) = vals(i);
    else
        x(i) = round(vals(i)/xsize);
        y(i) = rem(vals(i),xsize);
    end
end