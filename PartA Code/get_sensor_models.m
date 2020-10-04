clear, clc, close all

k = 1; % Figure number 

%% Import and Structure Data
filename = 'calibration.csv';
data = importdata(filename);

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

%% Remove sensor zero-readings
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

%% Visualise initial data
%{
figure(k)
k = k + 1;
scatter(data.range,data.sonar1);
title('sonar1')
ylabel('sensor reading')
xlabel('range')

figure(k)
k = k + 1;
scatter(data.range,data.sonar2);
title('sonar2')
ylabel('sensor reading')
xlabel('range')

figure(k)
k = k + 1;
scatter(data.range,data.raw_ir1);
title('ir1')
ylabel('sensor reading')
xlabel('range')

figure(k)
k = k + 1;
scatter(data.range,data.raw_ir2);
title('ir2')
ylabel('sensor reading')
xlabel('range')

figure(k)
k = k + 1;
scatter(data.range,data.raw_ir3);
title('ir3')
ylabel('sensor reading')
xlabel('range')
figure(k)

k = k + 1;
scatter(data.range,data.raw_ir4);
title('ir4')
ylabel('sensor reading')
xlabel('range')
%}

%% Remove Outliers & Find Coeffs
% Removal process is different for each sensor

% sonar1 outlier removal
tol = 0.04;
[clean.sonar1.range,clean.sonar1.signal,clean.sonar1.coeffs] =...
    rmv_outliers_poly(data.range,data.sonar1,1,tol);
sonar1.a = clean.sonar1.coeffs(1);
sonar1.b = clean.sonar1.coeffs(2);

% sonar2 outlier removal
tol = 0.04;
[clean.sonar2.range,clean.sonar2.signal,clean.sonar2.coeffs] =...
    rmv_outliers_poly(data.range,data.sonar2,1,tol);
sonar2.a = clean.sonar2.coeffs(1);
sonar2.b = clean.sonar2.coeffs(2);

% ir1 coeffs 
ir1_fit = fit(data.range,data.raw_ir1,'Rat21');
ir1.a = ir1_fit.p1; % 0.153;
ir1.b = ir1_fit.p2;
ir1.c = ir1_fit.p3;
ir1.d = ir1_fit.q1;

% ir2 coeffs (from cftool)
ir2_fit = fit(data.range,data.raw_ir2,'Rat21');
ir2.a = ir2_fit.p1; 
ir2.b = ir2_fit.p2;
ir2.c = ir2_fit.p3;
ir2.d = ir2_fit.q1;

% ir3 coeffs (from cftool)
ir3_fit = fit(data.range,data.raw_ir3,'Rat21');
ir3.a = ir3_fit.p1; 
ir3.b = ir3_fit.p2;
ir3.c = ir3_fit.p3;
ir3.d = ir3_fit.q1;

% ir4 coeffs (from cftool)
ir4_fit = fit(data.range,data.raw_ir4,'Rat32');
ir4.a = ir4_fit.p1; 
ir4.b = ir4_fit.p2;
ir4.c = ir4_fit.p3;
ir4.d = ir4_fit.p4;
ir4.e = ir4_fit.q1;
ir4.f = ir4_fit.q2;

% Combine into Clean Struct
sensor_coeffs.h.sonar1  = sonar1;
sensor_coeffs.h.sonar2  = sonar2;
sensor_coeffs.h.ir1     = ir1;
sensor_coeffs.h.ir2     = ir2;
sensor_coeffs.h.ir3     = ir3;
sensor_coeffs.h.ir4     = ir4;

%% Construct Sensor Models
% Breaks the data into segments (based on sensor reading) to calculate
% variance af a function of sensor reading.

x = data.range; % Rename for compactness

% Clearing outliers changed lengths of the sonar data, so have to be
% calculated seperately
sonar1_model = sonar1.a .* clean.sonar1.range + sonar1.b;
sonar2_model = sonar2.a .* clean.sonar2.range + sonar2.b;



ir1_model = (ir1.a.*x.^2+ir1.b.*x+ir1.c)./(x+ir1.d);
ir2_model = (ir2.a.*x.^2+ir2.b.*x+ir2.c)./(x+ir2.d);
ir3_model = (ir3.a.*x.^2+ir3.b.*x+ir3.c)./(x+ir3.d);
ir4_model = (ir4.a.*x.^3+ir4.b.*x.^2+ir4.c.*x+ir4.d)./(x.^2+ir4.e.*x+ir4.f);

%% Caluclate errors
sonar1_err = sonar1_model - clean.sonar1.signal;
sonar2_err = sonar2_model - clean.sonar2.signal;
ir1_err = ir1_model - data.raw_ir1;
ir2_err = ir2_model - data.raw_ir2;
ir3_err = ir3_model - data.raw_ir3;
ir4_err = ir4_model - data.raw_ir4;

% Combine into matrix for iteration: (sonars are different array lengths so
% are calculated seperately)

sonar1_errors(:,1)  = sonar1_err;
sonar1_errors(:,2)  = clean.sonar1.range;
sonar2_errors(:,1)  = sonar2_err;
sonar2_errors(:,2)  = clean.sonar2.range;
errors(:,1)         = ir1_err;
errors(:,2)         = ir2_err;
errors(:,3)         = ir3_err;
errors(:,4)         = ir4_err;
errors(:,5)         = data.range;

%% Calculate Variances

n = 50;        % Number of data points used for each segment

% ---------------------------------SONAR1-------------------------------
% Rearrange data points according to range
sonar1_errors = sortrows(sonar1_errors,2);
num_segs = floor(length(sonar1_errors)/n); % Last segment will include extra points
for j = 1:num_segs
    if j == num_segs
        % Special case, include the last remaining data points in final
        % segment
        segment = sonar1_errors((j-1)*n+1:length(sonar1_errors),1);
        range_segment = sonar1_errors((j-1)*n+1:length(sonar1_errors),2);
    else
        % Normal case
        segment = sonar1_errors((j-1)*n+1:j*n,1);
        range_segment = sonar1_errors((j-1)*n+1:j*n,2);
    end
    sonar1_variances(j,1) = var(segment);
    sonar1_variances(j,2) = mean(range_segment);
end

% ------------------------------SONAR2----------------------------------
% Rearrange data points according to range
sonar2_errors = sortrows(sonar2_errors,2);
num_segs = floor(length(sonar2_errors)/n); % Last segment will include extra points
for j = 1:num_segs
    if j == num_segs
        % Special case, include the last remaining data points in final
        % segment
        segment = sonar2_errors((j-1)*n+1:length(sonar2_errors),1);
        range_segment = sonar2_errors((j-1)*n+1:length(sonar2_errors),2);
    else
        % Normal case
        segment = sonar2_errors((j-1)*n+1:j*n,1);
        range_segment = sonar2_errors((j-1)*n+1:j*n,2);
    end
    sonar2_variances(j,1) = var(segment);
    sonar2_variances(j,2) = mean(range_segment);
end

% ----------------------------ALL IR SENSORS-----------------------------
for i = 1 : length(errors(1,:)) - 1
    % Rearrange data points according to range
    range_pos = length(errors(1,:));
    errors = sortrows(errors,range_pos);
    num_segs = floor(length(errors)/n); % Last segment will include extra points
    
    for j = 1:num_segs
        if j == num_segs
            % Special case, include the last remaining data points in final
            % segment
            segment = errors((j-1)*n+1:length(errors),i);
            range_segment = errors((j-1)*n+1:length(errors),range_pos);
        else
            % Normal case
            segment = errors((j-1)*n+1:j*n,i);
            range_segment = errors((j-1)*n+1:j*n,range_pos);
        end
        variances(j,i,1) = var(segment);
        variances(j,i,2) = mean(range_segment);
    end
end

%% Find Variance Coeffs

% sonar1
sonar1_var_fit = fit(sonar1_variances(:,2),sonar1_variances(:,1),'Poly2');
sensor_coeffs.v.sonar1.a = sonar1_var_fit.p1;
sensor_coeffs.v.sonar1.b = sonar1_var_fit.p2;
sensor_coeffs.v.sonar1.c = sonar1_var_fit.p3;

% sonar2
sonar2_var_fit = fit(sonar2_variances(:,2),sonar2_variances(:,1),'Poly2');
sensor_coeffs.v.sonar2.a = sonar2_var_fit.p1;
sensor_coeffs.v.sonar2.b = sonar2_var_fit.p2;
sensor_coeffs.v.sonar2.c = sonar2_var_fit.p3;

% ir1
ir1_var_fit = fit(variances(:,1,2),variances(:,1,1),'Rat21');
sensor_coeffs.v.ir1.a = ir1_var_fit.p1;
sensor_coeffs.v.ir1.b = ir1_var_fit.p2;
sensor_coeffs.v.ir1.c = ir1_var_fit.p3;
sensor_coeffs.v.ir1.d = ir1_var_fit.q1;

% ir2
ir2_var_fit = fit(variances(:,2,2),variances(:,2,1),'Rat21');
sensor_coeffs.v.ir2.a = ir2_var_fit.p1;
sensor_coeffs.v.ir2.b = ir2_var_fit.p2;
sensor_coeffs.v.ir2.c = ir2_var_fit.p3;
sensor_coeffs.v.ir2.d = ir2_var_fit.q1;

% ir3
ir3_var_fit = fit(variances(:,3,2),variances(:,3,1),'Poly4');
sensor_coeffs.v.ir3.a = ir3_var_fit.p1;
sensor_coeffs.v.ir3.b = ir3_var_fit.p2;
sensor_coeffs.v.ir3.c = ir3_var_fit.p3;
sensor_coeffs.v.ir3.d = ir3_var_fit.p4;
sensor_coeffs.v.ir3.e = ir3_var_fit.p5;

% ir4
ir4_var_fit = fit(variances(:,4,2),variances(:,4,1),'Rat21');
sensor_coeffs.v.ir4.a = ir4_var_fit.p1;
sensor_coeffs.v.ir4.b = ir4_var_fit.p2;
sensor_coeffs.v.ir4.c = ir4_var_fit.p3;
sensor_coeffs.v.ir4.d = ir4_var_fit.q1;

%% OLD SHIT INVERSE METHOD FOR VARIANCE (don't open pls) :(
%{
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
%}

%% SENSOR MODELS
% All coefficients and models used in this section were found using the
% get_sensor_models script, and using cftool to fit curves to the output
% sensor values and variance values
% --------------------------mean value models----------------------------
%sonar1
a = sensor_coeffs.h.sonar1.a;
b = sensor_coeffs.h.sonar1.b;
sonar1_h    = @(x) a.*x+b;
sonar1_dh   = @(x) a;

%sonar2
a = sensor_coeffs.h.sonar2.a;
b = sensor_coeffs.h.sonar2.b;
sonar2_h    = @(x) a.*x+b;
sonar2_dh   = @(x) a;

%ir1
a = sensor_coeffs.h.ir1.a;
b = sensor_coeffs.h.ir1.b;
c = sensor_coeffs.h.ir1.c;
d = sensor_coeffs.h.ir1.d;
ir1_h   = @(x) (a.*x.^2+b.*x+c)./(x+d);
ir1_dh  = @(x) ((x+d).*(2*a.*x+b)-(a.*x.^2+b.*x+c))./(x+d).^2;

%ir2
a = sensor_coeffs.h.ir2.a;
b = sensor_coeffs.h.ir2.b;
c = sensor_coeffs.h.ir2.c;
d = sensor_coeffs.h.ir2.d;
ir2_h   = @(x) (a.*x.^2+b.*x+c)./(x+d);
ir2_dh  = @(x) ((x+d).*(2*a.*x+b)-(a.*x.^2+b.*x+c))./(x+d).^2;

%ir3
a = sensor_coeffs.h.ir3.a;
b = sensor_coeffs.h.ir3.b;
c = sensor_coeffs.h.ir3.c;
d = sensor_coeffs.h.ir3.d;
ir3_h   = @(x) (a.*x.^2+b.*x+c)./(x+d);
ir3_dh  = @(x) ((x+d).*(2*a.*x+b)-(a.*x.^2+b.*x+c))./(x+d).^2;

%ir4
a = sensor_coeffs.h.ir4.a;
b = sensor_coeffs.h.ir4.b;
c = sensor_coeffs.h.ir4.c;
d = sensor_coeffs.h.ir4.d;
e = sensor_coeffs.h.ir4.e;
f = sensor_coeffs.h.ir4.f;
ir4_h   = @(x) (a*x.^3 + b*x.^2 + c*x + d)./(x.^2 + e*x + f);
ir4_dh  = @(x) ((x.^2+e.*x+f).*(3*a.*x.^2+2*b.*x+c)-(a.*x.^3+b.*x.^2+c.*x+d).*(2.*x+e))./(x.^2+e.*x+f).^2;


% ---------------------------variance models-----------------------------
% sonar1
a = sensor_coeffs.v.sonar1.a;
b = sensor_coeffs.v.sonar1.b;
c = sensor_coeffs.v.sonar1.c;
sonar1_V = @(x) a.*x.^2+b.*x+c;

%sonar2
a = sensor_coeffs.v.sonar2.a;
b = sensor_coeffs.v.sonar2.b;
c = sensor_coeffs.v.sonar2.c;
sonar2_V = @(x) a.*x.^2+b.*x+c;

%ir1
a = sensor_coeffs.v.ir1.a;
b = sensor_coeffs.v.ir1.b;
c = sensor_coeffs.v.ir1.c;
d = sensor_coeffs.v.ir1.d;
ir1_V = @(x) (a.*x.^2+b.*x+c)./(x+d);

%ir2
a = sensor_coeffs.v.ir2.a;
b = sensor_coeffs.v.ir2.b;
c = sensor_coeffs.v.ir2.c;
d = sensor_coeffs.v.ir2.d;
ir2_V = @(x) (a.*x.^2+b.*x+c)./(x+d);

%ir3
a = sensor_coeffs.v.ir3.a;
b = sensor_coeffs.v.ir3.b;
c = sensor_coeffs.v.ir3.c;
d = sensor_coeffs.v.ir3.d;
e = sensor_coeffs.v.ir3.e;
ir3_V = @(x) a.*x.^4+b.*x.^3+c.*x.^2+d.*x+e;

%ir4
a = sensor_coeffs.v.ir4.a;
b = sensor_coeffs.v.ir4.b;
c = sensor_coeffs.v.ir4.c;
d = sensor_coeffs.v.ir4.d;
ir4_V = @(x) (a.*x.^2+b.*x+c)./(x+d);

%% VISUALISE RESULTS
% Plotting Functions And Fit Lines
%
sonar1_fit = sonar1_h(data.range);
figure(k)
k = k+1;
subplot(2,3,1)
hold on
scatter(clean.sonar1.range,clean.sonar1.signal)
plot(data.range,sonar1_fit)
plot(data.range,sonar1_fit-sonar1_V(data.range).^0.25,'b')
plot(data.range,sonar1_fit+sonar1_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('sonar1')
hold off

sonar2_fit = sonar2_h(data.range);
subplot(2,3,2)
hold on
scatter(clean.sonar2.range,clean.sonar2.signal)
plot(data.range,sonar2_fit)
plot(data.range,sonar2_fit-sonar2_V(data.range).^0.25,'b')
plot(data.range,sonar2_fit+sonar2_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('sonar2')
hold off
%
ir1_fit = ir1_h(data.range);
subplot(2,3,3)
hold on
scatter(data.range,data.raw_ir1)
plot(data.range,ir1_fit)
plot(data.range,ir1_fit-ir1_V(data.range).^0.25,'b')
plot(data.range,ir1_fit+ir1_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('ir1')
hold off

ir2_fit = ir2_h(data.range);
subplot(2,3,4)
hold on
scatter(data.range,data.raw_ir2)
plot(data.range,ir2_fit)
plot(data.range,ir1_fit-ir2_V(data.range).^0.25,'b')
plot(data.range,ir2_fit+ir2_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('ir2')
hold off

ir3_fit = ir3_h(data.range);
subplot(2,3,5)
hold on
scatter(data.range,data.raw_ir3)
plot(data.range,ir3_fit)
plot(data.range,ir3_fit-ir3_V(data.range).^0.25,'b')
plot(data.range,ir3_fit+ir3_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('ir3')
hold off

ir4_fit = ir4_h(data.range);
subplot(2,3,6)
hold on
scatter(data.range,data.raw_ir4)
plot(data.range,ir4_fit)
plot(data.range,ir4_fit-ir4_V(data.range).^0.25,'b')
plot(data.range,ir4_fit+ir4_V(data.range).^0.25,'b')
xlabel('Range')
ylabel('Sensor Reading')
title('ir4')
hold off
%}

%% Plotting Variances

figure(k)
k = k + 1;
% Sonar1
subplot(2,3,1)
hold on
scatter(sonar1_variances(:,2),sonar1_variances(:,1));
plot(clean.sonar1.range,sonar1_V(clean.sonar1.range));
title('sonar1 variances')
xlabel('Range')
ylabel('Variance')
hold off

% Sonar2
x = clean.sonar2.range;
a = sensor_coeffs.v.sonar2.a;
b = sensor_coeffs.v.sonar2.b;
c = sensor_coeffs.v.sonar2.c;
sonar2_var_model = a.*x.^2+b.*x+c;
subplot(2,3,2)
hold on
scatter(sonar2_variances(:,2),sonar2_variances(:,1));
plot(clean.sonar2.range,sonar2_V(clean.sonar2.range));
title('sonar2 variances')
xlabel('Range')
ylabel('Variance')

% IR 1
subplot(2,3,3)
hold on
scatter(variances(:,1,2),variances(:,1,1));
plot(data.range,ir1_V(data.range));
title('ir1 variances')
xlabel('Range')
ylabel('Variance')

% IR 2
x = data.range;
a = sensor_coeffs.v.ir2.a;
b = sensor_coeffs.v.ir2.b;
c = sensor_coeffs.v.ir2.c;
d = sensor_coeffs.v.ir2.d;

ir2_var_model = (a.*x.^2+b.*x+c)/(x+d);
subplot(2,3,4)
hold on
scatter(variances(:,2,2),variances(:,2,1));
plot(data.range,ir2_V(data.range));
title('ir2 variances')
xlabel('Range')
ylabel('Variance')

% IR 3
subplot(2,3,5)
hold on
scatter(variances(:,3,2),variances(:,3,1));
plot(data.range,ir3_V(data.range));
title('ir3 variances')
xlabel('Range')
ylabel('Variance')

% IR 4
subplot(2,3,6)
hold on
scatter(variances(:,4,2),variances(:,4,1));
plot(data.range,ir4_V(data.range));
title('ir4 variances')
xlabel('Range')
ylabel('Variance')


