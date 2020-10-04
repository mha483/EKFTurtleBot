clear, clc

%% Import and Structure Data
filename = 'calibration.csv';
data = importdata(filename);

k = 1; % figure number

data.num = data.data(:,1);
data.time = data.data(:,2);
data.range = data.data(:,3);
data.vel_com = data.data(:,4);
data.raw_ir1 = data.data(:,5);
data.raw_ir2 = data.data(:,6);
data.raw_ir3 = data.data(:,7);
data.raw_ir4 = data.data(:,8);
data.sonar1 = data.data(:,9);
data.sonar2 = data.data(:,10);

% Put data into array for manipulation
dat_array(:,1) = data.range;
dat_array(:,2) = data.sonar1;
dat_array(:,3) = data.sonar2;
dat_array(:,4) = data.raw_ir1;
dat_array(:,5) = data.raw_ir2;
dat_array(:,6) = data.raw_ir3;
dat_array(:,7) = data.raw_ir4;

% Sort Data By Range
dat_array = sortrows(dat_array);

% Remove any rows where any sensor readings == 0
rowsToDel = max(dat_array==0,[],2);
dat_array(rowsToDel,:) = [];

% Move data back into struct
data.range      = dat_array(:,1);
data.sonar1     = dat_array(:,2);
data.sonar2     = dat_array(:,3);
data.raw_ir1    = dat_array(:,4);
data.raw_ir2    = dat_array(:,5);
data.raw_ir3    = dat_array(:,6);
data.raw_ir4    = dat_array(:,7);

% Visualise initial data
figure(k)
k = k + 1;
scatter(data.sonar1,data.range);
title('sonar1')
xlabel('sensor reading')
ylabel('range')
figure(k)
k = k + 1;
scatter(data.sonar2,data.range);
title('sonar2')
xlabel('sensor reading')
ylabel('range')
figure(k)
k = k + 1;
scatter(data.raw_ir1,data.range);
title('ir1')
xlabel('sensor reading')
ylabel('range')
figure(k)
k = k + 1;
scatter(data.raw_ir2,data.range);
title('ir2')
xlabel('sensor reading')
ylabel('range')
figure(k)
k = k + 1;
scatter(data.raw_ir3,data.range);
title('ir3')
xlabel('sensor reading')
ylabel('range')
figure(k)
k = k + 1;
scatter(data.raw_ir4,data.range);
title('ir4')
xlabel('sensor reading')
ylabel('range')


%% REMOVE OUTLIERS
% Equations found by removing outliers and fitting until 'good' fit is found
% Sonar 1
tol = 1;
[clean.sonar1.range,clean.sonar1.signal,clean.sonar1.coeffs] =...
    rmv_outliers_poly(data.range,data.sonar1,1);

[clean.sonar2.range,clean.sonar2.signal,clean.sonar2.coeffs] =...
    rmv_outliers_poly(data.range,data.sonar2,1);

[clean.ir1.range,clean.ir1.signal,clean.ir1.coeffs] =...
    rmv_outliers_pwr(data.range,data.raw_ir1,0.2);

[clean.ir2.range,clean.ir2.signal,clean.ir2.coeffs] =...
    rmv_outliers_pwr(data.range,data.raw_ir2,0.2);

[clean.ir3.range,clean.ir3.signal,clean.ir3.coeffs] =...
    rmv_outliers_pwr(data.range,data.raw_ir3,0.2);

[clean.ir4.range,clean.ir4.signal,clean.ir4.coeffs] =...
    rmv_outliers_pwr(data.range,data.raw_ir4,0.2);

%% FIND VARIANCES
% Breaks the data into segments (based on sensor reading) to calculate
% variance af a function of sensor reading.

n = 100; % Number of data points used for each segment

% Calculate Error Arrays
sonar1_range = (data.sonar1-clean.sonar1.coeffs(2))./clean.sonar1.coeffs(1);
sonar2_range = (data.sonar2-clean.sonar2.coeffs(2))./clean.sonar2.coeffs(1);
ir1_range = (data.raw_ir1./clean.ir1.coeffs(1)).^(1./clean.ir1.coeffs(2));
ir2_range = (data.raw_ir2./clean.ir2.coeffs(1)).^(1./clean.ir2.coeffs(2));
ir3_range = (data.raw_ir3./clean.ir3.coeffs(1)).^(1./clean.ir3.coeffs(2));
ir4_range = (data.raw_ir4./clean.ir4.coeffs(1)).^(1./clean.ir4.coeffs(2));

sonar1_err = sonar1_range - data.range;
sonar2_err = sonar2_range - data.range;
ir1_err = ir1_range - data.range;
ir2_err = ir2_range - data.range;
ir3_err = ir3_range - data.range;
ir4_err = ir4_range - data.range;

% Combine into matrix for iteration:
errors(:,1) = sonar1_err;
errors(:,2) = sonar2_err;
errors(:,3) = ir1_err;
errors(:,4) = ir2_err;
errors(:,5) = ir3_err;
errors(:,6) = ir4_err;
errors(:,7) = data.range;

for i = 1 : length(errors(1,:)) - 1
    % Rearrange data points according to sensor being considered
    errors = sortrows(errors,1);
    
    num_segs = floor(length(errors)/n); % Last segment will include extra points
    
    for j = 1:num_segs
        if j == num_segs
            % Special case, include the last remaining data points in final
            % segment
            segment = errors((j-1)*n+1:length(errors),i);
        else
            % Normal case
            segment = errors((j-1)*n+1:j*n,i);
        end
        variances(j,i) = var(segment);
    end
end

for i = 1:length(variances(1,:))
    figure(i+9)
    plot(1:length(variances(:,i)),variances(:,i));
end

%% VISUALISE RESULTS
% Non-inverted functions
%{
sonar1_fit = clean.sonar1.coeffs(1).*data.range + clean.sonar1.coeffs(2);
figure(2)
hold on
scatter(data.range,data.sonar1)
scatter(data.range,sonar1_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off

sonar2_fit = clean.sonar2.coeffs(1).*data.range + clean.sonar2.coeffs(2);
figure(3)
hold on
scatter(data.range,data.sonar2)
scatter(data.range,sonar2_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off

ir1_fit = clean.ir1.coeffs(1).*data.range.^clean.ir1.coeffs(2);
figure(4)
hold on
scatter(data.range,data.raw_ir1)
scatter(data.range,ir1_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off

ir2_fit = clean.ir2.coeffs(1).*data.range.^clean.ir2.coeffs(2);
figure(5)
hold on
scatter(data.range,data.raw_ir2)
scatter(data.range,ir2_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off

ir3_fit = clean.ir3.coeffs(1).*data.range.^clean.ir3.coeffs(2);
figure(6)
hold on
scatter(data.range,data.raw_ir3)
scatter(data.range,ir3_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off

ir4_fit = clean.ir4.coeffs(1).*data.range.^clean.ir4.coeffs(2);
figure(7)
hold on
scatter(data.range,data.raw_ir4)
scatter(data.range,ir4_fit)
xlabel('Range')
ylabel('Sensor Reading')
hold off
%}

% Inverted Functions
%{
sonar1_range = (data.sonar1-clean.sonar1.coeffs(2))./clean.sonar1.coeffs(1);
figure(2)
hold on
scatter(data.sonar1,data.range)
scatter(data.sonar1,sonar1_range)
title('Sonar1')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off

sonar2_range = (data.sonar2-clean.sonar2.coeffs(2))./clean.sonar2.coeffs(1);
figure(3)
hold on
scatter(data.sonar2,data.range)
scatter(data.sonar2,sonar2_range)
title('Sonar2')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off

ir1_range = (data.raw_ir1./clean.ir1.coeffs(1)).^(1./clean.ir1.coeffs(2));
figure(4)
hold on
scatter(data.raw_ir1,data.range)
scatter(data.raw_ir1,ir1_range)
title('Ir1')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off

ir2_range = (data.raw_ir2./clean.ir2.coeffs(1)).^(1./clean.ir2.coeffs(2));
figure(5)
hold on
scatter(data.raw_ir2,data.range)
scatter(data.raw_ir2,ir2_range)
title('Ir2')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off

ir3_range = (data.raw_ir3./clean.ir3.coeffs(1)).^(1./clean.ir3.coeffs(2));
figure(6)
hold on
scatter(data.raw_ir3,data.range)
scatter(data.raw_ir3,ir3_range)
title('Ir3')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off
%
ir4_range = (data.raw_ir4./clean.ir4.coeffs(1)).^(1./clean.ir4.coeffs(2));
figure(7)
hold on
scatter(data.raw_ir4,data.range)
scatter(data.raw_ir4,ir4_range)
title('Ir4')
xlabel('Sensor Reading')
ylabel('Range')
legend('True','Model')
hold off
%}
