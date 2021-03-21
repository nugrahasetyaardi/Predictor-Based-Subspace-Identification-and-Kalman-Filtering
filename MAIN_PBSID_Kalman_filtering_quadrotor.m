%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PBSID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contributor: Nugraha Setya Ardi - 10714522 - 952035
% Course: Estimation and Learning in Aerospace
% A/Y: 2020/2021


clc
clear
close all
load project_data
%% Initialization
Ts = 1/100;             % Sampling Time
t_in = t(1);            % Initial Sample
t_fin = t(end);         % Final Sample
t_interp = t;           % time interval
%% Getting dataset
u = setpoint_position_ned(:,2);
y1 = position_ned(:,2);
veloc = velocity_ned(:,2);

for i=1:length(t)
    q4 = q_body2ned(i,1);
    qvect = q_body2ned(i,2:end);
    q = [qvect q4];
    R = quat2dcm(q);
    a_ned0 = R*acceleration_body(i,:)';
    a_ned(i,:) = a_ned0';
end

y2 = a_ned(:,2);
figure
plot(t_interp,u)
figure
plot(t_interp,y1)
figure
plot(t_interp,y2)
ym = [y1 y2];
%% FFT and Filtering
y1f=fft_and_filter(y1,Ts,1,2,'position ned');
y2f=fft_and_filter(y2,Ts,1,2,'acceleration ned');
uf=fft_and_filter(u,Ts,3,2,'setpoint_position ned');
y1ff = [y1f y2f];
uf=u;
%% Define Identification and Validation Data Sets

% Create Entire Data Set
z = iddata(y1f,uf,Ts,'Name', 'Pitch Model', ...
    'InputName', 'setpoint_position_ned', 'OutputName', 'position_ned',...
    'InputUnit', 'm', 'OutputUnit', 'm'); advice(z);
% Create Entire Data Set
%% Plot subsets (2 Subsets)
samplesval=round(length(t)*2/3);                 % Samples per Subset
z1=z(1:samplesval);                           % Idenfitication Subset
z2=z(samplesval+1:end);                       % Validation Subset

figure;plot(z1,z2);legend('Idenfitication Subset','Validation Subset');

%% Get Data 
zz=get(z);
u_get=zz.InputData{1,1};    % Input
y_get=zz.OutputData{1,1};   % Filtered Output

%% --------------------- PBSID Implementation --------------------------%%

order1=5;       % Order of the model varx
f1=15;          % Future window varx
p1=f1;          % Past window varx

order2=5;       % Order of the model varmax
f2=15;          % Future window varmax
p2=f2;          % Past window varmax

u_1=uf(1:samplesval);                       % Input for identification
y_1=y1ff(1:samplesval,:);                   % Output for identification
u_val=uf(samplesval+1:end);                 % Input for validation
y_val=y1ff(samplesval+1:end,:);             % Output for validation
t_val = samplesval+1:length(t);             % Time for validation
t_val_plot = t(samplesval+1:length(t));     % Time for plot validation
y1f_val = y1f(samplesval+1:length(t));      % Position NED for validation
y2f_val = y2f(samplesval+1:length(t));      % Acceleration NED for validation

% figure
% plot(t(1:samplesval),uf(1:samplesval))
% hold on
% plot(t(samplesval+1:end),uf(samplesval+1:end))
% legend('Training dataset','Validation dataset')
% xlabel('time [s]')
% ylabel('setpoint position ned [m]')
% grid on
% 
% figure
% plot(t(1:samplesval),y1f(1:samplesval))
% hold on
% plot(t(samplesval+1:end),y1f(samplesval+1:end))
% legend('Training dataset','Validation dataset')
% xlabel('time [s]')
% ylabel('position ned [m]')
% grid on
% 
% figure
% plot(t(1:samplesval),y2f(1:samplesval))
% hold on
% plot(t(samplesval+1:end),y2f(samplesval+1:end))
% legend('Training dataset','Validation dataset')
% xlabel('time [s]')
% ylabel('acceleration ned [m/s^{-2}]')
% grid on


%% VARX
[S,X] = dordvarx(u_1,y_1,f1,p1,'tikh','gcv');
figure
subplot(2,1,1);
semilogy(S,'*');
title('Singular Values VARx');
x = dmodx(X,order1);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u_1,y_1,f1,p1);

SSvarx = ss(Ai,Bi,Ci,Di,1);
yvarx = lsim(SSvarx,u_val,t_val);
pos_varx = yvarx(:,1);
acc_varx = yvarx(:,2);

VAF_varx_pos = vaf(y1f_val,pos_varx);
VAF_varx_acc = vaf(y2f_val,acc_varx);

%% VARMAX
[S,X] = dordvarmax(u_1,y_1,f2,p2,'els',1e-6,'tikh','gcv');
subplot(2,1,2);
semilogy(S,'*');
title('Singular Values VARmax');
x = dmodx(X,order2);

[Av,Bv,Cv,Dv,Kv] = dx2abcdk(x,u_1,y_1,f2,p2);

SSvarmax = ss(Av,Bv,Cv,Dv,1);
yvarmax = lsim(SSvarmax,u_val,t_val);
pos_varmax = yvarmax(:,1);
acc_varmax = yvarmax(:,2);

VAF_varmax_pos = vaf(y1f_val,pos_varmax);
VAF_varmax_acc = vaf(y2f_val,acc_varmax);

%% Plotting
figure
plot(t_val,y1(samplesval+1:end),'--')
hold on
plot(t_val,y1f_val,'--');
hold on
plot(t_val,pos_varx);
hold on
plot(t_val,pos_varmax);
hold off
legend('measurement','smoothed measurement','VARX','VARMAX')

figure
plot(t_val,y2(samplesval+1:end),':')
hold on
plot(t_val,y2f_val,'--')
hold on
plot(t_val,acc_varx);
hold on
plot(t_val,acc_varmax);
hold off
legend('measurement','smoothed measurement','VARX','VARMAX')

%% Transfer function
[num,den] = ss2tf(Ai,Bi,Ci,Di);
TFvarxpos = tf(num(1,:),den,Ts);
TFvarxacc = tf(num(2,:),den,Ts);
[num,den] = ss2tf(Av,Bv,Cv,Dv);
TFvarmaxpos = tf(num(1,:),den,Ts);
TFvarmaxacc = tf(num(2,:),den,Ts);
z = iddata(y1f,uf,Ts);
zz = etfe(z);
z2 = iddata(y2f,uf,Ts);
zz2 = etfe(z2);

figure
bodemag(zz)
hold on
bodemag(TFvarxpos)
hold on
bodemag(TFvarmaxpos)
hold off
grid on
legend('smoothed measurement','VARX','VARMAX')
title('TF input to position')

figure
bodemag(zz2)
hold on
bodemag(TFvarxacc)
hold on
bodemag(TFvarmaxacc)
hold off
grid on
legend('smoothed measurement','VARX','VARMAX')
title('TF input to acceleration')

%% Plotting eigenvalues
figure
hold on
title('Poles of identified system (closed loop)')
[cx,cy] = pol2cart(linspace(0,2*pi),ones(1,100));
plot(cx,cy,'black');
plot(real(eig(Ai)),imag(eig(Ai)),'rx');
plot(real(eig(Av)),imag(eig(Av)),'bx');

axis([-1 1 -1 1]);
axis square
legend('STABBND','PBSID-varx','PBSID-varmax','Location','Best');
hold off

%% KALMAN DT-DT VARX

x0 = zeros(order1,1);   % initial state
P = eye(order1);        % initial state covariance
dt = 0.01;              % time step
tspan = t;              % span time

F = eye(order1)+ Ai*dt;
Q = 0.1*eye(order1);
H = Ci;
R = [1 0;0 1000]; 

for i=2:length(tspan)
    xhat(:,1) = x0;
    
    %prediction
    xmin = F*xhat(:,i-1) + Bi*u(i); 
    P = F*P*F'+Q;
    
    %correction
    K = P*H'*inv(H*P*H'+R);
    xhat(:,i) = xmin + K*(ym(i,:)'-H*xmin);
    P = (eye(order1)-K*H)*P;
    
    %output estimate
    yhat(:,i) = H*xhat(:,i) + Di*u(i);
    
    %lateral velocity ned
    velocity(i-1) = (yhat(1,i)-yhat(1,i-1))/dt;
end
velocity(length(tspan)) = velocity(length(tspan)-1);

%% KALMAN DT-DT VARMAX

x0 = zeros(order1,1);   % initial state
P = eye(order1);        % initial state covariance
dt = 0.01;              % time step
tspan = t;              % span time

F = eye(order2)+ Av*dt;
Q = 0.1*eye(order2);
H = Cv;
R = [1 0;0 1000]; 

for i=2:length(tspan)
    xhat(:,1) = x0;
    
    %prediction
    xmin = F*xhat(:,i-1) + Bv*u(i); 
    P = F*P*F'+Q;
    
    %correction
    K = P*H'*inv(H*P*H'+R);
    xhat(:,i) = xmin + K*(ym(i,:)'-H*xmin);
    P = (eye(order2)-K*H)*P;
    
    %output estimate
    yhat2(:,i) = H*xhat(:,i) + Dv*u(i);
    
    %lateral velocity ned
    velocity2(i-1) = (yhat(1,i)-yhat(1,i-1))/dt;
end
velocity2(length(tspan)) = velocity2(length(tspan)-1);

%% Plotting result

figure
plot(tspan,yhat(1,:));
hold on
plot(tspan,yhat2(1,:));
hold on
plot(tspan,y1);
hold on
plot(tspan,y1f);
hold off
legend('position - kalman filter VARX','position - kalman filter VARMAX','position - measurement','position - smoothed measurement')

figure
plot(tspan,y2,':');
hold on
plot(tspan,yhat(2,:));
hold on
plot(tspan,yhat2(2,:));
hold on
plot(tspan,y2f,'color',[0,0,0]);
hold off
legend('acceleration - measurement','acceleration - kalman filter VARX','acceleration - kalman filter VARMAX','acceleration - smoothed measurement')


figure
plot(tspan,velocity)
hold on
plot(tspan,velocity2)
hold on
plot(tspan,veloc,'color',[0,0,0])
hold off
legend('velocity - kalman filter VARX','velocity - kalman filter VARMAX','velocity dataset')

