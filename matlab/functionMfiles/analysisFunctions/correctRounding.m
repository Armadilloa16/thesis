function y = correctRounding(x,d)

if nargin == 2
    y = round((10^d)*x)/(10^d);
elseif nargin == 1
    y = round((1000)*x)/(1000);
else
    error('ERROR: !!! wrong number of inputs ')
end

end