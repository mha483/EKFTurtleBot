function [] = plotEKF(EKF_out)
% Plots the output of the EKF for the robotics assignment

time        = EKF_out.time;
x_post      = EKF_out.x_post;
v_post      = EKF_out.v_post;
x_prior     = EKF_out.x_prior;
v_prior     = EKF_out.v_prior;
sensor1     = EKF_out.sensor1;
sensor2     = EKF_out.sensor2;
sensor3     = EKF_out.sensor3;
w1          = EKF_out.w1;
w2          = EKF_out.w2;
w3          = EKF_out.w3;
K           = EKF_out.K;
range       = EKF_out.range;



figure(1)
subplot(4,1,1)
hold on
plot(time,x_post)
plot(time,x_post+v_post.^0.25,'r')
plot(time,x_post-v_post.^0.25,'r')
    
title('Posterior Values')
xlabel('time (s)')
ylabel('range (m)')
ylim([0,3.5])
hold off

subplot(4,1,2)

hold on
plot(time,v_post)
ylabel('Variance')


%{
hold on
plot(time,x_prior)
plot(time,x_prior+v_prior.^0.21,'r')
plot(time,x_prior-v_prior.^0.25,'r')
title('Prior Values')
xlabel('time (s)')
ylabel('range (m)')
ylim([0,3.5])
hold off
%}

subplot(4,1,3)
hold on
plot(time,w1)
plot(time,w2)
plot(time,w3)
legend(sensor1, sensor2, sensor3)
title('Sensor Weightings')
ylim([0,1])
hold off

subplot(4,1,4)
hold on
plot(time,K)
title('Kalman Gain')
ylim([0,1])


if ~isempty(range)
    figure(2)
    subplot(3,1,1)
    hold on
    plot(time,x_post)
    plot(time,range)
    legend('Filter', 'Error')
    title('Filter Comparison')
    ylabel('Range (m)')
    
    subplot(3,1,2)
    plot(time,x_post-range)
    ylabel('Error (m)')
    
    subplot(3,1,3)
    plot(time,v_post)
    ylabel('variance')
    xlabel('time')
end
end

