function cmap = divergingMap(n)
%DIVERGINGMAP Colorblind-conscious blue-white-orange diverging map.

arguments
    n (1,1) double {mustBeInteger,mustBePositive} = 257
end

blue = [0.12 0.31 0.50];
white = [0.96 0.96 0.94];
orange = [0.85 0.33 0.10];
half = floor(n / 2);
left = [linspace(blue(1), white(1), half + 1)', ...
    linspace(blue(2), white(2), half + 1)', ...
    linspace(blue(3), white(3), half + 1)'];
rightCount = n - half;
right = [linspace(white(1), orange(1), rightCount)', ...
    linspace(white(2), orange(2), rightCount)', ...
    linspace(white(3), orange(3), rightCount)'];
cmap = [left(1:end-1, :); right];
end
