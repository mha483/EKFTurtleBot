function [x,y] = clean_smallVals(x,y,tol)
% Remove entries <= 0, which will mess up power model.
while min(x) <= tol
    k = find(min(x) == x);
    y(k)   = [];
    x(k)    = [];
end

while min(y) <= tol
    k = find(min(y) == y);
    y(k)   = [];
    x(k)    = [];
end

end

