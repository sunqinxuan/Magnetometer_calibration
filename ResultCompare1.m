%{
�ֱ�֪ʶEKFУ����������ϣ��Ƚ϶��ߵ�Ч����

��дʱ�䣺
  2018.4.18
%}
%% format
clear;clc;
close all;

addpath('.\data')

data_original_filename = 'Flt1002_train.h5';
line_number = 1002.02; 

% data_info = h5info(data_original_filename);
data_line = h5read(data_original_filename,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

tt=readH5File(data_original_filename,'/tt',i1,i2);

flux_b_x=readH5File(data_original_filename,'/flux_b_x',i1,i2);
flux_b_y=readH5File(data_original_filename,'/flux_b_y',i1,i2);
flux_b_z=readH5File(data_original_filename,'/flux_b_z',i1,i2);
flux_b_t=readH5File(data_original_filename,'/flux_b_t',i1,i2);

w_x=readH5File(data_original_filename,'/roll_rate',i1,i2);
w_y=readH5File(data_original_filename,'/pitch_rate',i1,i2);
w_z=readH5File(data_original_filename,'/yaw_rate',i1,i2);


% flux_b_x = h5read(data_original_filename,'/flux_b_x');
% flux_b_x=flux_b_x(i1:i2,:);
% 
% flux_b_y = h5read(data_original_filename,'/flux_b_y');
% flux_b_y=flux_b_y(i1:i2,:);
% 
% flux_b_z = h5read(data_original_filename,'/flux_b_z');
% flux_b_z=flux_b_z(i1:i2,:);
% 
% flux_b_t = h5read(data_original_filename,'/flux_b_t');
% flux_b_t=flux_b_t(i1:i2,:);


%% load the data
% [gyro_still] = xlsread('./stableB.xlsx');                                   % ��ֹʱ�����ǵ����������
% figure
% plot(gyro_still(1:1000,4:6))

% [data_test] = xlsread('huawei_x1.xlsx');
% mag_test = data_test(:,1:3);
mag_test = [flux_b_x,flux_b_y,flux_b_z]*0.01; % convert from nT to mG;

% [data_test]  = xlsread ('���ݼ�\fastwalking_swing_circle.xlsx');     
% mag_test = data_test(:,1:3)/10;
testlength      = size( mag_test, 1 );

% [M]          = xlsread ('���ݼ�\fastwalking_swing_circle.xlsx');                          % �����ʼ����
% data.b_p     = M(:,1:3)/10;                                                    % �����Ʋ���ֵ����λΪmG
% data.w       = M(:,4:6)/1800*pi;           
% [M]          = xlsread ('huawei_x1.xlsx');                          % �����ʼ����
% data.b_p     = M(:,1:3);                                                    % �����Ʋ���ֵ����λΪmG
% data.w       = M(:,4:6);                                                    % �����ǲ���ֵ����λrad/s
% data.dt      = 0.01;                                                        % ����ʱ��������λs
data.b_p     = [flux_b_x,flux_b_y,flux_b_z]*0.01; % convert from nT to mG;
data.w       = [w_x,w_y,w_z]*pi/180.0;            % convert from deg/s to rad/s;
data.dt      = 0.1;
data.mrw     = 0.5;                                                         % �ų�b��random walks����λΪmG������������table I���õ�ֵ
data.wrw     = 0.1/180*pi;                                                  % ���ٶ�w��random walks,��λΪdegree/��s^1/2������������table I���õ�ֵ                                     % ��������
data.m       = size( data.b_p, 1 );
data.phi     = 50;                                                         % �����й�ʽ17�еĲ��������ڱ�ʾ�����ǵĿɿ��̶ȡ�
data.P       = [500*eye(3) zeros(3,6)  zeros(3);
                zeros(6,3) 1e-4*eye(6) zeros(6,3);
                zeros(3)   zeros(3,6)  500*eye(3)];                         % �����������������table I���õ�ֵ
W = eye(3);
V = [0 0 0]';
%% calibrated gyroscope measurement at time k-1��equ.16
% w_bias     = mean(gyro_still(1:1000,4:6))/1800*pi;
% for j = 1: data.m
%     w_noise(j,:)     = w_bias; 
% end
w            = data.w;
sz_w         = size(w);
w_noise      = data.wrw * randn(sz_w); 
w_cal        = data.w - w_noise;

%% Select correction interval
% 5s movement data for calibration
% start       = 500;
% range       = start+1:1:start+1+1000;
range       = 1:1:testlength-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF Calibration
k = 1;
for i = range 
    %% make the datain
    datain.h_p1       = data.b_p(i+1,:)';
    datain.h_p0       = data.b_p(i,:)';
    datain.dt         = data.dt;
    datain.phi        = data.phi;
    datain.mrw        = data.mrw;
    datain.wrw        = data.wrw;   
    datain.w_cal      = w_cal(i,:);
    if k == 1 
    % Setting the initial value
        datain.W      = eye(3);
        datain.V      = [0 0 0]';
        datain.P      = data.P;
        datain.h_p0   = datain.W ^(-1) * (datain.h_p0 - datain.V);
    else
    % make the value equal to the value at last time 
        datain.W      = dataout.W;
        datain.V      = dataout.V;
        datain.P      = dataout.P;
        datain.h_p0   = dataout.B;
    end
    
    %% EKF
    [dataout] = fun_MagCal_EKF(datain);
    
    %% load the output
    B_cal_EKF(k,:)    = dataout.B;
    W_cal_EKF(:,:,k)  = [dataout.W(1,1) dataout.W(1,2) dataout.W(1,3);
                         dataout.W(1,2) dataout.W(2,2) dataout.W(2,3);
                         dataout.W(1,3) dataout.W(2,3) dataout.W(3,3)];
    V_cal_EKF(k,:)    = dataout.V;
    P_cal_EKF(k,:)    = [dataout.P(1,1) dataout.P(4,4) dataout.P(5,5) dataout.P(10,10)];
    h_pre_EKF(k,:)    = dataout.h_p_pre;
    Kk(:,:,k)    = dataout.Kk;

    k = k + 1;
end
%% output
EKF_W = W_cal_EKF(:,:,end);
EKF_V =V_cal_EKF(end,:)';

figure 
plot(V_cal_EKF(:,1))
% ��ʾ���
fprintf( '��ž��� EKF_W:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n',EKF_W);
fprintf( 'Ӳ��ʸ�� EKF_V:\n  %g %g %g\n',EKF_V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ellipsoid Fitting
h = 1;
for i = range 
    dataElli(h,:) = data.b_p(i,:);
    h = h + 1; 
end
[Elli_Winv, Elli_V, Elli_B, Elli_E] = f_Mag_Ellipsoid_Fit(dataElli,10);
% ��ʾ���
fprintf( '��ž��� Elli_Winv:\n  %g %g %g\n  %g %g %g\n  %g %g %g\n',Elli_Winv);
fprintf( 'Ӳ��ʸ�� Elli_V:\n  %g %g %g\n',Elli_V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare
for i= 1:testlength
    normdata(i,:)      = norm(mag_test(i,:));

    B_EKF_WV(i,:)      = EKF_W^(-1)*(mag_test(i,:)' - EKF_V);
    normdata_EKF(i,:)  = norm(B_EKF_WV(i,:));

    B_Elli_WV(i,:)     = Elli_Winv * (mag_test(i,:)' - Elli_V);
    normdata_Elli(i,:) = norm(B_Elli_WV(i,:));

end


% figure
% plot3(h_pre_EKF(:,1),h_pre_EKF(:,2),h_pre_EKF(:,3),'.')
% hold on
% plot3(data.b_p(range,1),data.b_p(range,2),data.b_p(range,3),'.')
% hold on 
% % plot3(B_cal_EKF(:,1),B_cal_EKF(:,2),B_cal_EKF(:,3),'.')
% axis equal
% grid on 



%% comppute the magnetic field strength
Bfield_EKF  = mean(normdata_EKF);
Bfield_Elli = mean(normdata_Elli);
Bfield_raw      = mean(normdata);
standard = 50 * ones(testlength ,1);

quality1_EKF  = std(normdata_EKF) / Bfield_EKF
quality1_Elli = std(normdata_Elli) / Bfield_Elli

difference_raw = normdata-standard;
difference_EKF = normdata_EKF-standard;
difference_Elli = normdata_Elli-standard;
figure
plot(difference_raw,'r');hold on 
plot(difference_EKF,'b');hold on 
plot(difference_Elli,'g');hold on 
legend('ԭʼ����','EKF������','Elli������')
xlabel('ʱ��(unit:0.01s)')
ylabel('����(unit:mG)')



quality2_raw  = mean(normdata-standard)
quality2_EKF  = mean(normdata_EKF-standard)
quality2_Elli = mean(normdata_Elli-standard)

figure
grid on
% plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'.')
% hold on
 plot3(B_EKF_WV(:,1),B_EKF_WV(:,2),B_EKF_WV(:,3),'.')
% hold on
[x, y, z] = ellipsoid(0,0,0,Bfield_EKF,Bfield_EKF,Bfield_EKF);
mesh(x, y, z, 'FaceAlpha', 0.1, 'EdgeColor', 'k');
hold on
[~, X0, ~] = EllipsoidFitting(mag_test(:,1),mag_test(:,2),mag_test(:,3), 1, 1);
[~, X_EKF, ~] = EllipsoidFitting( B_EKF_WV(:,1) , B_EKF_WV(:,2) , B_EKF_WV(:,3), 1, 2);
[~, X_Elli, ~] = EllipsoidFitting( B_Elli_WV(:,1),B_Elli_WV(:,2),B_Elli_WV(:,3), 1, 3);
legend('����ų�����','ԭʼ����','EKF����','Elli����');
% legend('����ų�����','ԭʼ����','EKF����');

xlabel('X (unit:mG)')
ylabel('Y (unit:mG)')
zlabel('Z (unit:mG)')
plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'r.');hold on
plot3(data.b_p(range,1),data.b_p(range,2),data.b_p(range,3),'ys');hold on
axis vis3d;
axis equal;
grid on
% title('Ч��չʾ')
% axis([-80 80 -80 80 -80 80 ])

% axis equal;
% grid on
% title('Elli')
% [x_Elli, y_Elli, z_Elli] = ellip% figure
% plot3(data.b_p(:,1),data.b_p(:,2),data.b_p(:,3),'.')
% hold on
% plot3(B_Elli_WV(:,1),B_Elli_WV(:,2),B_Elli_WV(:,3),'.')
% axis vis3d;
% ellipsoid(0,0,0,Bfield_Elli,Bfield_Elli,Bfield_Elli);
% mesh(x_Elli, y_Elli, z_Elli, 'FaceAlpha', 0.1, 'EdgeColor', 'r');





