function [range,signal,coeffs] = rmv_outliers_poly(range,signal,n,tol)
% Removes outliers from a data set and returns clean data set
%   range   = true range values
%   signal  = signal values
%   n       = degree of polynomial used to approximate curve
%   tol     = Outlier rejection ends when error is below tolerance

num = 10;   % Number of points deleted each iteration


for i = 1:9999
    
    a = polyfit(range,signal,n);
    line_fit = 0;
    % Construct the polynomial approximation curve
    for h = 1:length(a)
        line_fit = line_fit + range.^(length(a)-h) * a(h);
    end
    % Error in current poly. approximation for each point
    error = abs(line_fit-signal);
    
    % Visualise progress with each iteration
    %{
    plot(range,line_fit);
    hold on
    scatter(range,signal);
    hold off
    pause(0.1)
    %}
    
    % Find and delete *num* data points with greatest error (outliers)
    for j = 1:num
        k = find(max(error) == error);
        error(k)    = [];
        signal(k)   = [];
        range(k)    = [];
    end
 
    if max(error) < tol
        break
    end
end
coeffs = a;
end

