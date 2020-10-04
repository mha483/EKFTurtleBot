function [x,y,coeffs] = rmv_outliers_pwr(x,y,tol)
% Removes outliers from a data set and returns clean data set
% Outlier rejection ends when error is below tolerance
%   range   = true range values
%   signal  = signal values

num = 10;   % Number of points deleted each iteration

% Remove entries <= 0, which will mess up power model.
while min(x) <= 0
    k = find(min(x) == x);
    y(k)   = [];
    x(k)    = [];
end

while min(y) <= 0
    k = find(min(y) == y);
    y(k)   = [];
    x(k)    = [];
end



for i = 1:9999
    
    a = fit(x,y,'power1');
    line_fit = a.a*x.^(a.b);
    % Error in current approximation for each point
    error = abs(line_fit-y);
    
    % Visualise progress with each iteration
    %{
    scatter(x,line_fit);
    hold on
    scatter(x,y);
    hold off
    %}
    
    % Find and delete *num* data points with greatest error (outliers)
    for j = 1:num
        k = find(max(error) == error);
        error(k)    = [];
        y(k)   = [];
        x(k)    = [];
    end
 
    if max(error) < tol
        break
    end
end
coeffs = [a.a,a.b];
end

