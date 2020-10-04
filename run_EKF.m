% Kalman filter to estimate robot position
% For ENMT481
% Author:       Max Harrison
% Last Edited:  01/09/20

clear, clc

%% IMPORT DATA
filename = 'training2.csv';
data = importdata(filename);
if length(data.data(1,:)) == 10
    i = 1; % Range is known
elseif length(data.data(1,:)) == 9
    i = 0;
else
    error('Bad Input');
end

filename = 'sensorModelCoeffs.mat';
sensor_coeffs = importdata(filename);

filename = 'motionModelCoeffs.mat';
motion_coeffs = importdata(filename);

time        = data.data(:,2);
vel_com     = data.data(:,3+i);    
ir1         = data.data(:,4+i);
ir2         = data.data(:,5+i);
ir3         = data.data(:,6+i); 
ir4         = data.data(:,7+i);
sonar1      = data.data(:,8+i);
sonar2      = data.data(:,9+i);

if i
    range = data.data(:,3);
else
    range = [];
end


%% DEFINITIONS

mean_x_post = 0.1; % Initial Pos
var_x_post  = 1; % Initial Certainty

% Initialise arrays

EKF_out.time        = time(1:end-1,1);
EKF_out.x_prior     = zeros(length(time)-1,1);
EKF_out.v_prior     = zeros(length(time)-1,1);
EKF_out.x_post      = zeros(length(time)-1,1);
EKF_out.v_post      = zeros(length(time)-1,1);
EKF_out.sonar1_x    = zeros(length(time)-1,1);
EKF_out.sonar2_x    = zeros(length(time)-1,1);
EKF_out.ir1_x       = zeros(length(time)-1,1);
EKF_out.ir2_x       = zeros(length(time)-1,1);
EKF_out.ir3_x       = zeros(length(time)-1,1);
EKF_out.ir4_x       = zeros(length(time)-1,1);
EKF_out.sonar1_var  = zeros(length(time)-1,1);
EKF_out.sonar2_var  = zeros(length(time)-1,1);
EKF_out.ir1_var     = zeros(length(time)-1,1);
EKF_out.ir2_var     = zeros(length(time)-1,1);
EKF_out.ir3_var     = zeros(length(time)-1,1);
EKF_out.ir4_var     = zeros(length(time)-1,1);
EKF_out.w_sonar1    = zeros(length(time)-1,1);
EKF_out.w_ir1       = zeros(length(time)-1,1);
EKF_out.w_ir4       = zeros(length(time)-1,1);
EKF_out.fused_x     = zeros(length(time)-1,1);
EKF_out.fused_var   = zeros(length(time)-1,1);
EKF_out.K           = zeros(length(time)-1,1);

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
k = 1;     % Scalable gain to reduce trust of sonar1
% sonar1
a = sensor_coeffs.v.sonar1.a;
b = sensor_coeffs.v.sonar1.b;
c = sensor_coeffs.v.sonar1.c;
sonar1_V = @(x) (a.*x.^2+b.*x+c).*k;

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

%% MOTION MODELS
% All coefficients and models used in this section were found using the
% get_vehicle_model script, which generates a function to estimate velocity
% based on the commanded velocity
% --------------------------mean value model----------------------------
a = motion_coeffs.h.a;
b = motion_coeffs.h.b;
est_v_mean = @(x) a.*x + b;
% ---------------------------variance model-----------------------------
a = motion_coeffs.v.a;
b = motion_coeffs.v.b;
c = motion_coeffs.v.c;
est_v_var = @(x) a.*x.^2 + b.*x + c;

%% EXTENDED KALMAN FILTER

for i = 1:length(time)-1
    
    dt = time(i+1) - time(i);
    
    % Estimate actual velocity (motion model)
    est_vel_mean = est_v_mean(vel_com(i));
    % Get estimate's variance (var_W)
    est_vel_var = est_v_var(vel_com(i)) * dt^2;
    
    % Calculate prior
    mean_x_prior  = mean_x_post + est_vel_mean * dt;
    var_x_prior   = var_x_post + est_vel_var;
    
    % Linearise sensor models based on prior estimate
    % Invert sensor models
    % Use sensor models to estimate positions (x1_hat, x2_hat, x3_hat)
    sonar1_x = (sonar1(i) - sonar1_h(mean_x_prior))/sonar1_dh(mean_x_prior) + mean_x_prior;
    sonar2_x = (sonar2(i) - sonar2_h(mean_x_prior))/sonar2_dh(mean_x_prior) + mean_x_prior;
    ir1_x = (ir1(i) - ir1_h(mean_x_prior))/ir1_dh(mean_x_prior) + mean_x_prior;
    ir2_x = (ir2(i) - ir2_h(mean_x_prior))/ir2_dh(mean_x_prior) + mean_x_prior;
    ir3_x = (ir3(i) - ir3_h(mean_x_prior))/ir3_dh(mean_x_prior) + mean_x_prior;
    ir4_x = (ir4(i) - ir4_h(mean_x_prior))/ir4_dh(mean_x_prior) + mean_x_prior;
    % Get sensor variances too
    sonar1_var = sonar1_V(mean_x_prior)   / sonar1_dh(mean_x_prior)^2;
    sonar2_var = sonar2_V(mean_x_prior)   / sonar2_dh(mean_x_prior)^2;
    ir1_var   = ir1_V(mean_x_prior)       / ir1_dh(mean_x_prior)^2;
    ir2_var   = ir2_V(mean_x_prior)       / ir2_dh(mean_x_prior)^2;
    ir3_var   = ir3_V(mean_x_prior)       / ir3_dh(mean_x_prior)^2;
    ir4_var   = ir4_V(mean_x_prior)       / ir4_dh(mean_x_prior)^2;
    
    % Choose 3 sensors
    sensor1 = 'sonar1';
    sensor2 = 'ir2';
    sensor3 = 'ir4';
    
    h1 = sonar1_x;
    h2 = ir2_x;
    h3 = ir4_x;
    
    v1 = sonar1_var;
    v2 = ir2_var;
    v3 = ir4_var;
    
    %
    % Find outliers in sonar sensors and ignore them
    if abs(h1-mean_x_prior) > 0.1 && h1==sonar1_x || h1==sonar2_x
        v1 = 999; % = very high number, reading is essentially ignored
    end
    
    if abs(h2-mean_x_prior) > 0.1 && h2==sonar1_x || h2==sonar2_x
        v2 = 999; % = very high number, reading is essentially ignored
    end
    
    if abs(h3-mean_x_prior) > 0.1 && h3==sonar1_x || h3==sonar2_x
        v3 = 999; % = very high number, reading is essentially ignored
    end
    %}
    
    % Fuse sensor readings
    w1          = (1/v1)/(1/v1+1/v2+1/v3);
    w2          = (1/v2)/(1/v1+1/v2+1/v3);
    w3          = (1/v3)/(1/v1+1/v2+1/v3);
    
    fused_x     = w1*h1 + w2*h2 + w3*h3;
    fused_var   = 1/(1/v1 + 1/v2 + 1/v3);
    
    % Calculate Kalman Gain
    K = 1/(fused_var) / (1/fused_var + 1/var_x_prior);
    
    % Calculate Posterior
    mean_x_post = mean_x_prior + K * (fused_x - mean_x_prior);
    var_x_post = (1 - K) * var_x_prior;
    
    % Store Values
    EKF_out.x_prior(i)      = mean_x_prior;
    EKF_out.v_prior(i)      = var_x_prior;
    EKF_out.x_post(i)       = mean_x_post;
    EKF_out.v_post(i)       = var_x_post;
    EKF_out.sonar1_x(i)     = sonar1_x;
    EKF_out.sonar2_x(i)     = sonar2_x;
    EKF_out.ir1_x(i)        = ir1_x;
    EKF_out.ir2_x(i)        = ir2_x;
    EKF_out.ir3_x(i)        = ir3_x;
    EKF_out.ir4_x(i)        = ir4_x;
    EKF_out.sonar1_var(i)   = sonar1_var;
    EKF_out.sonar2_var(i)   = sonar2_var;
    EKF_out.ir1_var(i)      = ir1_var;
    EKF_out.ir2_var(i)      = ir2_var;
    EKF_out.ir3_var(i)      = ir3_var;
    EKF_out.ir4_var(i)      = ir4_var;
    EKF_out.w1(i)           = w1;
    EKF_out.w2(i)           = w2;
    EKF_out.w3(i)           = w3;
    EKF_out.fused_x(i)      = fused_x;
    EKF_out.fused_var(i)    = fused_var;
    EKF_out.K(i)            = K;
    
    
end

if k ~= 1
    fprintf('CAUTION: sonar1 variance was scaled by k = %.0f\n', k)
end

EKF_out.sensor1      = sensor1;
EKF_out.sensor2      = sensor2;
EKF_out.sensor3      = sensor3;
EKF_out.range        = range(2:end);

plotEKF(EKF_out);



