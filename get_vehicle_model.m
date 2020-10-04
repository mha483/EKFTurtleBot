clear, clc, close all

k = 1; % Figure number 

%% Import and Structure Data
filename = 'training2.csv';
data = importdata(filename);

time    = data.data(:,2);
range   = data.data(:,3);
vel_com = data.data(1:end-1,4); % Remove last entry to match estimated vel

%% Estimate Actual Velocity

vel = zeros(length(time)-1,1);

for i = 1:length(time)-1
    dt = time(i+1)  - time(i);
    dx = range(i+1) - range(i);
    vel(i) = dx/dt;
end

%% Find Motion Model Mean

a = polyfit(vel_com,vel,1); % Linear fit

% Put into struct
motion_coeffs.h.a = a(1);
motion_coeffs.h.b = a(2);

% Estimated velocities based on the motion model
vel_est = a(1).*vel_com + a(2);

%% Find Motion Model Variance

n = 30; % Number of data points used for each segment

errors(:,1) = vel_est - vel;
errors(:,2) = vel_com;
errors = sortrows(errors,2); % Sort errors by commanded velocity (the independent variable)
num_segs = floor(length(errors)/n); % Last segment will include up to n-1 extra points

% Caluculate variances for each segment
for j = 1:num_segs
    if j == num_segs
        % Special case, include the last remaining data points in final
        % segment
        segment = errors((j-1)*n+1:length(errors),1);
        range_segment = errors((j-1)*n+1:length(errors),2);
    else
        % Normal case
        segment = errors((j-1)*n+1:j*n,1);
        range_segment = errors((j-1)*n+1:j*n,2);
    end
    variances(j,1) = var(segment);
    variances(j,2) = mean(range_segment);
end

b = polyfit(variances(:,2),variances(:,1),2);
% Put into struct
motion_coeffs.v.a = b(1);
motion_coeffs.v.b = b(2);
motion_coeffs.v.c = b(3);
% create fit curve
var_est = b(1).*variances(:,2).^2 + b(2).*variances(:,2)+b(3);



%% Visualise Results
figure(1)
hold on
title('Motion Model Fitting')
xlabel('Commanded Velocity')
ylabel('Actual Velocity')
scatter(vel_com,vel)
plot(vel_com,vel_est)
hold off

figure(2)
hold on
title('Motion Model Variance Fitting')
xlabel('Commanded Velocity')
ylabel('Variance')
scatter(variances(:,2),variances(:,1))
plot(variances(:,2),var_est)
hold off

fprintf('Model mean coeffs are a=%.4d and b=%.4d\n',a(1),a(2))
fprintf('Model var coeffs are a = %.4d and b = %.4d and c = %.4d\n',b(1),b(2),b(3))

%% Comparison with Other Dataset
%{
% Storing previous values for comparison
vel_com_old = vel_com;
vel_est_old = vel_est;
var_est_old = var_est;
variances_old = variances;

errors = [];

% Import and Structure Data
filename = 'training1.csv';
data = importdata(filename);

time    = data.data(:,2);
range   = data.data(:,3);
vel_com = data.data(1:end-1,4); % Remove last entry to match estimated vel

% Estimate Actual Velocity

vel = zeros(length(time)-1,1);

for i = 1:length(time)-1
    dt = time(i+1)  - time(i);
    dx = range(i+1) - range(i);
    vel(i) = dx/dt;
end

% Find Motion Model Mean

a = polyfit(vel_com,vel,1); % Linear fit

% Estimated velocities based on the motion model
vel_est = a(1).*vel_com + a(2);

% Find Motion Model Variance

n = 30; % Number of data points used for each segment

errors(:,1) = vel_est - vel;
errors(:,2) = vel_com;
errors = sortrows(errors,2); % Sort errors by commanded velocity (the independent variable)
num_segs = floor(length(errors)/n); % Last segment will include up to n-1 extra points

% Caluculate variances for each segment
for j = 1:num_segs
    if j == num_segs
        % Special case, include the last remaining data points in final
        % segment
        segment = errors((j-1)*n+1:length(errors),1);
        range_segment = errors((j-1)*n+1:length(errors),2);
    else
        % Normal case
        segment = errors((j-1)*n+1:j*n,1);
        range_segment = errors((j-1)*n+1:j*n,2);
    end
    variances(j,1) = var(segment);
    variances(j,2) = mean(range_segment);
end

b = polyfit(variances(:,2),variances(:,1),2);
var_est = b(1).*variances(:,2).^2 + b(2).*variances(:,2)+b(3);

% Visualise Results
figure(3)
hold on
title('Comparison With Alternative Training Data')
xlabel('Commanded Velocity')
ylabel('Actual Velocity')
scatter(vel_com,vel)
plot(vel_com,vel_est)
plot(vel_com_old,vel_est_old)
legend('Training1 Data','Training1 Fit','Training2 Fit')
hold off

figure(4)
hold on
title('Motion Model Variance Fitting')
xlabel('Commanded Velocity')
ylabel('Variance')
scatter(variances(:,2),variances(:,1))
plot(variances(:,2),var_est)
plot(variances_old(:,2),var_est_old)
legend('Training1 Variances','Training1 Variance Function','Training2 Variance Function')
hold off
%}




