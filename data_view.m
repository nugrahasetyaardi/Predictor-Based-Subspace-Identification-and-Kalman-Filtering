clear,clc
close all

load project_data

% Plot position setpoint, position (NED frame), acceleration (body frame),
% velocity (NED frame)

figure(1)
subplot(411)
plot(t,setpoint_position_ned)
grid
ylabel('Position setpoint [m]')

subplot(412)
plot(t,position_ned)
grid
ylabel('Position [m]')

subplot(413)
plot(t,acceleration_body)
grid
ylabel('Acceleration [m/s^2]')

subplot(414)
plot(t,velocity_ned)
grid
xlabel(' time [s]')
ylabel('Velocity [m/s]')

% Plot quaternion corresponding to rotation from body to NED (convention:
% q=[q4 q1 q2 q3]

figure(2)
subplot(411)
plot(t,q_body2ned(:,1))
grid

subplot(412)
plot(t,q_body2ned(:,2))
grid

subplot(413)
plot(t,q_body2ned(:,3))
grid

subplot(414)
plot(t,q_body2ned(:,4))
grid
xlabel(' time [s]')

