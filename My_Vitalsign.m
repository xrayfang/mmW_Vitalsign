%% 采用硬件平台 ：AWR1243EVM+DCA1000EVM
%% 时间：2024年07月28日
%% 功能：单人呼吸心跳原始数据采集与生命体征信号处理与提取
%% 采集环境：办公桌前，胸部正对采集设备0.58m距离
%% 算法流程：预处理+MTI+反正切提取相位+相位解缠+相位差分+滑动平均去噪+带通滤波器+提取估计包络归一化心跳波
%% ========================================================================
clc;clear;close all;
%% =========================================================================
%% 雷达参数设置
c = physconst('lightspeed');
f0 = 77E9;
lambda = c/f0;
ADCStartTime = 6E-6;
IdleTime = 10E-6;
RampEndTime = 63.14E-6;
K = 63.343E12;               % 调频率 (Hz/sec)
fStart = f0+K*ADCStartTime;  % 起始频率
fS = 9121e3;                 % ADC采样频率 (sps)
Ts = 1/fS;                   % ADC采样间隔
Nchirp = 2400;                % 一个帧的chirp数
Nadc = 512;                  % 一个chirp的adc采样点数
B = Nadc*Ts*K;
rangeRes = c/(2*B);
Nrx = 4;                     % 接收天线Rx数量
Ntx = 1;                     % 发射天线Tx数量
Nant = Nrx*Ntx;
%% 读取数据部分
% 读取数据
filename = 'adc_data.bin';
adcData = readDCA1000(filename);
%% 重排回波数据
rawDataTx = zeros(Nrx*Ntx,Nadc,Nchirp); 
% format: [Nrx*Ntx,Nadc,Nchirp]
% [Rx0Tx0 Rx0Tx1 | Rx1Tx0 Rx1Tx1 | Rx2Tx0 Rx2Tx1 | Rx3Tx0 Rx3Tx1]
channal = 0;    %初始化通道数
tic
for ii = 1:Nrx
    tempData = squeeze(adcData(ii,:));               
    tempData = reshape(tempData,Nadc,[]);
    for jj = 1:Ntx
        channal = channal + 1; % 通道数递增1 最大通道数channalMax = Nrx*Ntx
        rawDataTx(channal,:,:)  = tempData(:,jj:Ntx:end);
    end
end
toc
% [Rx0Tx0 Rx1Tx0 Rx2Tx0 Rx3Tx0 | Rx0Tx1 Rx1Tx1 Rx2Tx1 Rx3Tx1]
rawDataRx = zeros(Nrx*Ntx,Nadc,Nchirp);
for jj = 1:Ntx
    rawDataRx(1+Nrx*(jj-1):Nrx*jj,:,:) = rawDataTx(jj:Ntx:end,:,:);
end
% --------------------------------------------------------------------------------
% 虚拟 阵元设置 根据Rx和Tx的位置关系设置
% [Rx3Tx0 Rx2Tx0 Rx1Tx0 Rx0Tx0 | Rx3Tx1 Rx2Tx1 Rx1Tx1 Rx0Tx1]
rawDataVx = zeros(Nrx*Ntx,Nadc,Nchirp);
for jj = 1:Ntx
    rawDataVx(1+Nrx*(jj-1):Nrx*jj,:,:) = rawDataTx(flip(jj:Ntx:end),:,:);
end
rawData = squeeze(rawDataVx(1,:,:));
%% =============经过以上处理获得了最为关键的等待信号算法处理的数据：rawData（512*1024矩阵）==============
%% 距离维FFT
% 时延误差校正
tI = 5.6248e-10; % Instrument delay for range calibration (corresponds to a 8.44cm range offset)
t = (0:Nadc-1)*Ts;
rangeBiasCorrector = exp(-1j*2*pi*K*tI*t).';
rawData = rawData.*rangeBiasCorrector;
% 距离FFT
nFFTtime = 2400;
rangeData = fft(rawData,nFFTtime).';
% 作图
rangeAxis = rangeRes*(0:nFFTtime-1)*(Nadc/nFFTtime);
[X,Y] = meshgrid(rangeAxis,(1:Nchirp));  
figure;
mesh(X,Y,abs(rangeData));
xlabel('距离(m)');
ylabel('脉冲chrip数');
zlabel('幅度');
title('距离维-1DFFT结果');

%% 消除静态干扰
% 非必须
% buffer = zeros(2400,2400);
% buffer(2:end,:) =  rangeData(1:end-1,:);
% rangeData = rangeData - buffer;
% 
% figure;
% mesh(X,Y,abs(rangeData));
% xlabel('距离(m)');
% ylabel('脉冲chrip数');
% zlabel('幅度');
% title('距离维-1DFFT结果(MTI)');

%% 提取相位
angleData = angle(rangeData);

detaR = rangeRes*(Nadc/nFFTtime);
% Range-bin tracking 找出能量最大的点，即人体的位置  
rangeAbs = abs(rangeData);

rangeSum = sum(rangeAbs);
% 限定检测距离
minRange = 0.3;
maxRange = 1.5;
minIndex = floor(minRange/detaR);
maxIndex = ceil(maxRange/detaR);
rangeSum(1:minIndex) = 0;
rangeSum(maxIndex:end) = 0;
[~,index] = max(rangeSum);

%% 取出能量最大点的相位  extract phase from selected range bin
angleTarget = angleData(:,index);%1024个chrip的某列range bin的相位
% 提取相位信号（原始）
figure;
plot(1:Nchirp,angleTarget);
xlabel('时间/点数（N）：对应每个chrip');
ylabel('相位');
title('未展开相位信号');

phi=angleTarget;
%% 进行相位解缠
% unwrap函数：
phi=unwrap(phi);
figure;
plot(phi);
xlabel('点数（N）：对应每个chrip');
ylabel('相位（rad）');
title('解缠后的相位');
angle_fft_last = phi;

%% phase difference 相位差分
%通过减去连续的相位值，对展开的相位执行相位差运算，
angle_fft_last2=zeros(1,Nchirp);
for i = 2:Nchirp
    angle_fft_last2(i) = angle_fft_last(i) - angle_fft_last(i-1);
end 

figure;
plot(angle_fft_last2);
xlabel('点数（N）：对应每个chrip');
ylabel('相位（rad）');
title('相位差分后信号');

%% 脉冲噪声去除：滑动平均滤波（选取0.25s(50ms*5)的滑动窗口，窗口长度为5）
% 去除由于测试环境引起的脉冲噪声
phi=smoothdata(angle_fft_last2,'movmean',5);
figure;
plot(phi);
title('滑动平均滤波相位信号');
%对相位差分信号作FFT 
N1=length(phi);
FS=20;                            % 帧周期50ms --> Fs = 1/0.05 = 20Hz
FFT = abs(fft(phi));              %--FFT取模，幅度
f=(0:N1-1)*(FS/N1);             %其中每点的频率
%傅里叶变换结果对称
figure;
plot(f(1:N1/8),FFT(1:N1/8)) %取前一部分放大观察
xlabel('频率（f/Hz）');
ylabel('幅度');
title('相位信号FFT  ');
%%  IIR带通滤波 Bandpass Filter 0.1-0.5hz，输出呼吸信号
fs =20; %呼吸心跳信号的采样率
%设计IIR，4阶巴特沃斯带通滤波器
load('coe3.mat', 'breath_pass');
breath_data = filter(breath_pass,phi); 

figure;
plot(breath_data);
xlabel('时间(s)');
ylabel('幅度');
title('呼吸时域波形');
%% FFT-Peak
breath = abs(fft(phi));                                     

figure;
plot(f(1:130),breath(1:130));   %取前一部分放大观察
xlabel('频率（f/Hz）');
ylabel('幅度');
title('呼吸信号FFT  ');

[~,breath_index] = max(breath); %谱峰最大值搜索

breath_count =(fs*(breath_index-1)/Nchirp)*60;        %呼吸频率解算
%% IIR带通滤波 Bandpass Filter 0.8-2hz，输出心跳的数据
%  设计IIR，8阶巴特沃斯带通滤波器
load('coe4.mat', 'heart_pass');

heart_data = filter(heart_pass,phi); 
 
figure;
plot(heart_data);
xlabel('时间(s)');
ylabel('幅度');
title('心跳时域波形');

N1=length(heart_data);
f=(0:N1-1)*(fs/N1);       %其中每点的频率
heart=abs(fft(heart_data));

figure;
plot(f(1:200),heart(1:200));
xlabel('频率（f/Hz）');
ylabel('幅度');
title('心跳信号FFT');

[~,heart_index] = max(heart);
% 判断是否有人
%
heart_count =(fs*(heart_index-1)/Nchirp)*60;%心跳频率解算
%% 1024帧，51.2s
% 如果数据长度够长，则雷达会51.2s对呼吸数据和心跳数据进行一次刷新
%以便实现更为精确的检测。
disp(['呼吸：',num2str(breath_count),'  心跳：',num2str(heart_count)])

%% END 
