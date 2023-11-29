try
    radar1.stop();
    radar1.close();
    radar2.stop();
    radar2.close();
    radar3.stop();
    radar3.close();
    radar4.stop();
    radar4.close();
    radar5.stop();
    radar5.close();
    radar6.stop();
    radar6.close();
    radar7.stop();
    radar7.close();
    radar8.stop();
    radar8.close();
catch me
end

clear;
close all;

%% Port Number
COMP4 = 'COM18';   % port number
COMP5 = 'COM19';   % port number
COMP6 = 'COM20';   % port number
COMP7 = 'COM21';   % port number
COMP8 = 'COM14';   % port number
COMP9 = 'COM15';   % port number
COMP0 = 'COM16';   % port number
COMP1 = 'COM17';   % port number


%% Radar Setup
FPS = 20;
CollectionH = 1/6;
SaveOnOffFlag = 1;
TxPower = 3;
TxCenterFreq = 4;
FrameStart = 0.4; % meters.
FrameStop = 7; % meters.
CorrectMeter = 2;

i1 = 0;
DataCursor1 = 1;
ResetCounter1 = 0;
dataType = 'rf';
% chip setting
PPS = 9*100; % 400
DACmin = 849+75;
DACmax = 1200-75;
Iterations = 4;
SamplingFreq = 23.328; % GHz
one_img_fr = 64;

%bandpass setting
if(TxCenterFreq==3)
    CutOff_Low = 6.0;
    CutOff_High = 8.5;
    Bandwidth = (CutOff_High - CutOff_Low); % GHz
    CenterFreq = 7.290; % GHz
elseif(TxCenterFreq==4)    
    CutOff_Low = 7.25;
    CutOff_High = 10.20;
    Bandwidth = (CutOff_High - CutOff_Low); % GHz
    CenterFreq = 8.748; % GHz
end
N    = 40;       % Order
Fc1  = CutOff_Low/(SamplingFreq/2);      % First Cutoff Frequency
Fc2  = CutOff_High/(SamplingFreq/2);     % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);
% Calculate the coefficients using the FIR1 function.
bpf  = fir1(N, [Fc1 Fc2], 'bandpass', win, flag);

%% Load Libraries
addpath('matlab');
addpath('include');
addpath('lib');

Lib = ModuleConnector.Library;

%% Create X4radar object
radar1 = BasicRadarClassX4(COMP4,FPS,dataType);

radar1.open();
radar1.init();

% Configure X4 chip.
radar1.radarInstance.x4driver_set_pulsesperstep(PPS);
radar1.radarInstance.x4driver_set_dac_min(DACmin);
radar1.radarInstance.x4driver_set_dac_max(DACmax);
radar1.radarInstance.x4driver_set_iterations(Iterations);
radar1.radarInstance.x4driver_set_tx_power(TxPower);
radar1.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar1.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart1, frameStop1] = radar1.radarInstance.x4driver_get_frame_area();

radar1.start();

%chip2
radar2 = BasicRadarClassX4(COMP5,FPS,dataType);

radar2.open();
radar2.init();

% Configure X4 chip.
radar2.radarInstance.x4driver_set_pulsesperstep(PPS);
radar2.radarInstance.x4driver_set_dac_min(DACmin);
radar2.radarInstance.x4driver_set_dac_max(DACmax);
radar2.radarInstance.x4driver_set_iterations(Iterations);
radar2.radarInstance.x4driver_set_tx_power(TxPower);
radar2.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar2.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart2, frameStop2] = radar2.radarInstance.x4driver_get_frame_area();

radar2.start();

%chip3
radar3 = BasicRadarClassX4(COMP6,FPS,dataType);

radar3.open();
radar3.init();

% Configure X4 chip.
radar3.radarInstance.x4driver_set_pulsesperstep(PPS);
radar3.radarInstance.x4driver_set_dac_min(DACmin);
radar3.radarInstance.x4driver_set_dac_max(DACmax);
radar3.radarInstance.x4driver_set_iterations(Iterations);
radar3.radarInstance.x4driver_set_tx_power(TxPower);
radar3.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar3.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart3, frameStop3] = radar3.radarInstance.x4driver_get_frame_area();

radar3.start();


%chip4
radar4 = BasicRadarClassX4(COMP7,FPS,dataType);

radar4.open();
radar4.init();

% Configure X4 chip.
radar4.radarInstance.x4driver_set_pulsesperstep(PPS);
radar4.radarInstance.x4driver_set_dac_min(DACmin);
radar4.radarInstance.x4driver_set_dac_max(DACmax);
radar4.radarInstance.x4driver_set_iterations(Iterations);
radar4.radarInstance.x4driver_set_tx_power(TxPower);
radar4.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar4.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart4, frameStop4] = radar4.radarInstance.x4driver_get_frame_area();

radar4.start();

%% Create X4radar object
radar5 = BasicRadarClassX4(COMP8,FPS,dataType);

radar5.open();
radar5.init();

% Configure X4 chip.
radar5.radarInstance.x4driver_set_pulsesperstep(PPS);
radar5.radarInstance.x4driver_set_dac_min(DACmin);
radar5.radarInstance.x4driver_set_dac_max(DACmax);
radar5.radarInstance.x4driver_set_iterations(Iterations);
radar5.radarInstance.x4driver_set_tx_power(TxPower);
radar5.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar5.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart5, frameStop5] = radar5.radarInstance.x4driver_get_frame_area();

radar5.start();

%chip2
radar6 = BasicRadarClassX4(COMP9,FPS,dataType);

radar6.open();
radar6.init();

% Configure X4 chip.
radar6.radarInstance.x4driver_set_pulsesperstep(PPS);
radar6.radarInstance.x4driver_set_dac_min(DACmin);
radar6.radarInstance.x4driver_set_dac_max(DACmax);
radar6.radarInstance.x4driver_set_iterations(Iterations);
radar6.radarInstance.x4driver_set_tx_power(TxPower);
radar6.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar6.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart6, frameStop6] = radar6.radarInstance.x4driver_get_frame_area();

radar6.start();

%chip3
radar7 = BasicRadarClassX4(COMP0,FPS,dataType);

radar7.open();
radar7.init();

% Configure X4 chip.
radar7.radarInstance.x4driver_set_pulsesperstep(PPS);
radar7.radarInstance.x4driver_set_dac_min(DACmin);
radar7.radarInstance.x4driver_set_dac_max(DACmax);
radar7.radarInstance.x4driver_set_iterations(Iterations);
radar7.radarInstance.x4driver_set_tx_power(TxPower);
radar7.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar7.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart7, frameStop7] = radar7.radarInstance.x4driver_get_frame_area();

radar7.start();


%chip4
radar8 = BasicRadarClassX4(COMP1,FPS,dataType);

radar8.open();
radar8.init();

% Configure X4 chip.
radar8.radarInstance.x4driver_set_pulsesperstep(PPS);
radar8.radarInstance.x4driver_set_dac_min(DACmin);
radar8.radarInstance.x4driver_set_dac_max(DACmax);
radar8.radarInstance.x4driver_set_iterations(Iterations);
radar8.radarInstance.x4driver_set_tx_power(TxPower);
radar8.radarInstance.x4driver_set_tx_center_frequency(TxCenterFreq);
% Configure frame area
radar8.radarInstance.x4driver_set_frame_area(FrameStart,FrameStop);
% Read back actual set frame area
[frameStart8, frameStop8] = radar8.radarInstance.x4driver_get_frame_area();

radar8.start();
tstart = tic;  % time log start

%% figure option
f = figure;
f.Position = [50 50 1300 700];

subplot(3,3,1);
p1 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th1 = title('');

subplot(3,3,2);
p2 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th2 = title('');

subplot(3,3,3);
p3 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th3 = title('');

subplot(3,3,4);
p4 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th4 = title('');

subplot(3,3,5);
p5 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th5 = title('');

subplot(3,3,6);
p6 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th6 = title('');

subplot(3,3,7);
p7 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th7 = title('');

subplot(3,3,8);
p8 = plot(0,0);
xlabel('Distance [cm]');
ylabel('Amplitude');
ylim([0 0.05]);
xlim([0 200]);

grid on;
th8 = title('');

rl = 40;
rl2 = rl/2;
XYZ_beacons=[0 rl rl;0 rl 0;0 0 rl;0 0 0;rl rl rl;rl rl 0;rl 0 0;rl 0 rl];
initial_guess=[rl2 rl2 rl2];

f2 = figure;
p9 = plot3(0,0,0,'g');
p9.LineWidth = 1;
xlim([0 rl]);
ylim([0 rl]);
zlim([0 rl]);
grid on
th9 = title('');

%% User Variable Setup
Alpha = 0.95;
Timefor1Write = 3.2;  % Time for one hand-writing [sec]
TotalDataNum = 10;  % Total number of Data in one file 

OneDataLength = FPS*Timefor1Write;            % Frame number of one data
%MaxDataLength = OneDataLength*TotalDataNum;   % Total Frame number
MaxDataLength = 1;

position = [0,0,0];
record_length = 20;
record_posit = zeros(record_length,3);

%% Radar Measurment Loop
try
    while(1)
        numPackets1 = radar1.bufferSize();         % Peek message data float

        if (numPackets1 > 0)
            %% Get Raw Frame 
            i1 = i1+1;
            % Get frame (uses read_message_data_float)
            [frame1, ctr] = radar1.GetFrameNormalized();
            frame_Raw1 = frame1;
            [frame2, ctr] = radar2.GetFrameNormalized();
            frame_Raw2 = frame2;
            [frame3, ctr] = radar3.GetFrameNormalized();
            frame_Raw3 = frame3;
            [frame4, ctr] = radar4.GetFrameNormalized();
            frame_Raw4 = frame4;
            [frame5, ctr] = radar5.GetFrameNormalized();
            frame_Raw5 = frame5;
            [frame6, ctr] = radar6.GetFrameNormalized();
            frame_Raw6 = frame6;
            [frame7, ctr] = radar7.GetFrameNormalized();
            frame_Raw7 = frame7;
            [frame8, ctr] = radar8.GetFrameNormalized();
            frame_Raw8 = frame8;
            
            
            
            %% Data Processing
            
            if ( (DataCursor1==1)&&(ResetCounter1==0) ) % initialize only at the first iteration.
                DataLength1 = length(frame1);
                DataLengthCorrected1 = round(DataLength1*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori1 = zeros(MaxDataLength,DataLengthCorrected1);  % Original frame
                RawData1 = zeros(MaxDataLength,DataLengthCorrected1);  % RawData
                Hilbert_History1 = zeros(MaxDataLength,DataLengthCorrected1); % Hilbert Transform History
                
                ClockArray1 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter1 = zeros(1,DataLengthCorrected1);  % clutter
                B_LowPass1 = zeros(MaxDataLength,DataLengthCorrected1); % After Clutter Removed      
                
                %radar2
                DataLength2 = length(frame2);
                DataLengthCorrected2 = round(DataLength2*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori2 = zeros(MaxDataLength,DataLengthCorrected2);  % Original frame
                RawData2 = zeros(MaxDataLength,DataLengthCorrected2);  % RawData
                Hilbert_History2 = zeros(MaxDataLength,DataLengthCorrected2); % Hilbert Transform History
                
                ClockArray2 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter2 = zeros(1,DataLengthCorrected2);  % clutter
                B_LowPass2 = zeros(MaxDataLength,DataLengthCorrected2); % After Clutter Removed   

                %radar3
                DataLength3 = length(frame3);
                DataLengthCorrected3 = round(DataLength3*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori3 = zeros(MaxDataLength,DataLengthCorrected3);  % Original frame
                RawData3 = zeros(MaxDataLength,DataLengthCorrected3);  % RawData
                Hilbert_History3 = zeros(MaxDataLength,DataLengthCorrected3); % Hilbert Transform History
                
                ClockArray3 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter3 = zeros(1,DataLengthCorrected3);  % clutter
                B_LowPass3 = zeros(MaxDataLength,DataLengthCorrected3); % After Clutter Removed 

                %radar4
                DataLength4 = length(frame4);
                DataLengthCorrected4 = round(DataLength4*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori4 = zeros(MaxDataLength,DataLengthCorrected4);  % Original frame
                RawData4 = zeros(MaxDataLength,DataLengthCorrected4);  % RawData
                Hilbert_History4 = zeros(MaxDataLength,DataLengthCorrected4); % Hilbert Transform History
                
                ClockArray4 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter4 = zeros(1,DataLengthCorrected4);  % clutter
                B_LowPass4 = zeros(MaxDataLength,DataLengthCorrected4); % After Clutter Removed 
                   
                %radar5
                DataLength5 = length(frame5);
                DataLengthCorrected5 = round(DataLength5*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori5 = zeros(MaxDataLength,DataLengthCorrected5);  % Original frame
                RawData5 = zeros(MaxDataLength,DataLengthCorrected5);  % RawData
                Hilbert_History5 = zeros(MaxDataLength,DataLengthCorrected5); % Hilbert Transform History
                
                ClockArray5 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter5 = zeros(1,DataLengthCorrected5);  % clutter
                B_LowPass5 = zeros(MaxDataLength,DataLengthCorrected5); % After Clutter Removed      
                
                %radar6
                DataLength6 = length(frame6);
                DataLengthCorrected6 = round(DataLength6*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori6 = zeros(MaxDataLength,DataLengthCorrected6);  % Original frame
                RawData6 = zeros(MaxDataLength,DataLengthCorrected6);  % RawData
                Hilbert_History6 = zeros(MaxDataLength,DataLengthCorrected6); % Hilbert Transform History
                
                ClockArray6 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter6 = zeros(1,DataLengthCorrected6);  % clutter
                B_LowPass6 = zeros(MaxDataLength,DataLengthCorrected6); % After Clutter Removed   

                %radar7
                DataLength7 = length(frame7);
                DataLengthCorrected7 = round(DataLength7*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori7 = zeros(MaxDataLength,DataLengthCorrected7);  % Original frame
                RawData7 = zeros(MaxDataLength,DataLengthCorrected7);  % RawData
                Hilbert_History7 = zeros(MaxDataLength,DataLengthCorrected7); % Hilbert Transform History
                
                ClockArray7 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter7 = zeros(1,DataLengthCorrected7);  % clutter
                B_LowPass7 = zeros(MaxDataLength,DataLengthCorrected7); % After Clutter Removed 

                %radar8
                DataLength8 = length(frame8);
                DataLengthCorrected8 = round(DataLength8*CorrectMeter/(FrameStop - FrameStart)); 
                RawData_Ori8 = zeros(MaxDataLength,DataLengthCorrected8);  % Original frame
                RawData8 = zeros(MaxDataLength,DataLengthCorrected8);  % RawData
                Hilbert_History8 = zeros(MaxDataLength,DataLengthCorrected8); % Hilbert Transform History
                
                ClockArray8 = zeros(MaxDataLength,6);   % Clock Array
                LowPass_Clutter8 = zeros(1,DataLengthCorrected8);  % clutter
                B_LowPass8 = zeros(MaxDataLength,DataLengthCorrected8); % After Clutter Removed 
            end
            
            %radar1
            RawData_Ori1(1,:) = frame1(1:DataLengthCorrected1);
            RawData1(1,:) = RawData_Ori1(1,:);
            RawData1(1,:) = RawData1(1,:) - mean(RawData1(1,:));      % DC Removal
            RawData1(1,:) = conv(RawData1(1,:),bpf,'same');                       % Bandpass Filtering
            
            LowPass_Clutter1 = Alpha .* LowPass_Clutter1 + (1-Alpha) .* RawData1(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass1(1,:) = RawData1(1,:) - LowPass_Clutter1;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp1 =  abs(hilbert(B_LowPass1(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History1(DataCursor1,:) = Hilbert_temp1;                                        % Make Hilbert History 

            %radar2
            RawData_Ori2(1,:) = frame2(1:DataLengthCorrected2);
            RawData2(1,:) = RawData_Ori2(1,:);
            RawData2(1,:) = RawData2(1,:) - mean(RawData2(1,:));      % DC Removal
            RawData2(1,:) = conv(RawData2(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter2 = Alpha .* LowPass_Clutter2 + (1-Alpha) .* RawData2(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass2(1,:) = RawData2(1,:) - LowPass_Clutter2;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp2 =  abs(hilbert(B_LowPass2(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History2(DataCursor1,:) = Hilbert_temp2;                                        % Make Hilbert History 

            %radar3
            RawData_Ori3(1,:) = frame3(1:DataLengthCorrected3);
            RawData3(1,:) = RawData_Ori3(1,:);
            RawData3(1,:) = RawData3(1,:) - mean(RawData3(1,:));      % DC Removal
            RawData3(1,:) = conv(RawData3(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter3 = Alpha .* LowPass_Clutter3 + (1-Alpha) .* RawData3(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass3(1,:) = RawData3(1,:) - LowPass_Clutter3;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp3 =  abs(hilbert(B_LowPass3(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History3(DataCursor1,:) = Hilbert_temp3;                                        % Make Hilbert History 

            %radar4
            RawData_Ori4(1,:) = frame4(1:DataLengthCorrected4);
            RawData4(1,:) = RawData_Ori4(1,:);
            RawData4(1,:) = RawData4(1,:) - mean(RawData4(1,:));      % DC Removal
            RawData4(1,:) = conv(RawData4(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter4 = Alpha .* LowPass_Clutter4 + (1-Alpha) .* RawData4(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass4(1,:) = RawData4(1,:) - LowPass_Clutter4;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp4 =  abs(hilbert(B_LowPass4(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History4(DataCursor1,:) = Hilbert_temp4;                                        % Make Hilbert History 

            %radar5
            RawData_Ori5(1,:) = frame5(1:DataLengthCorrected5);
            RawData5(1,:) = RawData_Ori5(1,:);
            RawData5(1,:) = RawData5(1,:) - mean(RawData5(1,:));      % DC Removal
            RawData5(1,:) = conv(RawData5(1,:),bpf,'same');                       % Bandpass Filtering
            
            LowPass_Clutter5 = Alpha .* LowPass_Clutter5 + (1-Alpha) .* RawData5(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass5(1,:) = RawData5(1,:) - LowPass_Clutter5;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp5 =  abs(hilbert(B_LowPass5(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History1(DataCursor1,:) = Hilbert_temp1;                                        % Make Hilbert History 

            %radar6
            RawData_Ori6(1,:) = frame6(1:DataLengthCorrected6);
            RawData6(1,:) = RawData_Ori6(1,:);
            RawData6(1,:) = RawData6(1,:) - mean(RawData6(1,:));      % DC Removal
            RawData6(1,:) = conv(RawData6(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter6 = Alpha .* LowPass_Clutter6 + (1-Alpha) .* RawData6(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass6(1,:) = RawData6(1,:) - LowPass_Clutter6;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp6 =  abs(hilbert(B_LowPass6(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History2(DataCursor1,:) = Hilbert_temp2;                                        % Make Hilbert History 

            %radar7
            RawData_Ori7(1,:) = frame7(1:DataLengthCorrected7);
            RawData7(1,:) = RawData_Ori7(1,:);
            RawData7(1,:) = RawData7(1,:) - mean(RawData7(1,:));      % DC Removal
            RawData7(1,:) = conv(RawData7(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter7 = Alpha .* LowPass_Clutter7 + (1-Alpha) .* RawData7(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass7(1,:) = RawData7(1,:) - LowPass_Clutter7;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp7 =  abs(hilbert(B_LowPass7(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History3(DataCursor1,:) = Hilbert_temp3;                                        % Make Hilbert History 

            %radar8
            RawData_Ori8(1,:) = frame8(1:DataLengthCorrected8);
            RawData8(1,:) = RawData_Ori8(1,:);
            RawData8(1,:) = RawData8(1,:) - mean(RawData8(1,:));      % DC Removal
            RawData8(1,:) = conv(RawData8(1,:),bpf,'same');                       % Bandpass Filtering
            
            
            LowPass_Clutter8 = Alpha .* LowPass_Clutter8 + (1-Alpha) .* RawData8(1,:);    % Background Subtraction (Clutter Update)            
            B_LowPass8(1,:) = RawData8(1,:) - LowPass_Clutter8;                 % Background Subtraction (Clutter Removal)
            
            Hilbert_temp8 =  abs(hilbert(B_LowPass8(1,:)));                              % Hilbert Transform & abs function
            %Hilbert_History4(DataCursor1,:) = Hilbert_temp4;                                        % Make Hilbert History 

            FDataCursor1 = DataCursor1;
            DataCursor1 = DataCursor1 + 1;
                       
            
        end
        

        %% Update Plot
        set(p1, 'XData', (1:DataLengthCorrected1)/1.5625, 'YData', Hilbert_temp1)
        set(p2, 'XData', (1:DataLengthCorrected2)/1.5625, 'YData', Hilbert_temp2)
        set(p3, 'XData', (1:DataLengthCorrected3)/1.5625, 'YData', Hilbert_temp3)
        set(p4, 'XData', (1:DataLengthCorrected4)/1.5625, 'YData', Hilbert_temp4)
        set(p5, 'XData', (1:DataLengthCorrected5)/1.5625, 'YData', Hilbert_temp5)
        set(p6, 'XData', (1:DataLengthCorrected6)/1.5625, 'YData', Hilbert_temp6)
        set(p7, 'XData', (1:DataLengthCorrected7)/1.5625, 'YData', Hilbert_temp7)
        set(p8, 'XData', (1:DataLengthCorrected8)/1.5625, 'YData', Hilbert_temp8)
        th1.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th2.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th3.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th4.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th5.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th6.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th7.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        th8.String = strcat("RawData Envelope - Frame No : ",num2str(FDataCursor1));
        [amp1, arg1] = max(Hilbert_temp1);
        [amp2, arg2] = max(Hilbert_temp2);
        [amp3, arg3] = max(Hilbert_temp3);
        [amp4, arg4] = max(Hilbert_temp4);
        [amp5, arg5] = max(Hilbert_temp5);
        [amp6, arg6] = max(Hilbert_temp6);
        [amp7, arg7] = max(Hilbert_temp7);
        [amp8, arg8] = max(Hilbert_temp8);

        measured_dist = [arg1; arg2; arg3; arg4; arg5; arg6; arg7; arg8];
        measured_amp = [amp1; amp2; amp3; amp4; amp5; amp6; amp7; amp8];
        threshold = 0.0045;
        alpha = 1;
        % survive_dist=[];
        % survive_anch=[];
        % survive_posit=[];
        % 
        index = measured_amp>threshold;
        survive_dist = measured_dist(index);
        survive_anch = XYZ_beacons(index,:);

        if sum(index, "all")<4
            position = position;
        else
            old_position = position;
            dist_fun=@(pos) sqrt(sum((bsxfun(@minus,survive_anch,pos(:)')).^2,2));
            cost_function=@(pos) sum((dist_fun(pos)-survive_dist).^2);
            position=fminsearch(cost_function,initial_guess);
            %position=fminsearchbnd(cost_function,initial_guess, [50 50 50], [0 0 0]);

            if position(1)<0 || position(2)<0 || position(3)<0 || position(1)>rl || position(2)>rl || position(3)>rl
                position = old_position;
            end

            % comb_index = nchoosek(1:length(survive_dist),4);
            % length_index = size(comb_index);
            % survive_posit = zeros(length_index(1),3);
            % 
            % for i=1:length_index(1)
            %     dist_fun=@(pos) sqrt(sum((bsxfun(@minus,survive_anch(comb_index(i,:), :),pos(:)')).^2,2));
            %     cost_function=@(pos) sum((dist_fun(pos)-survive_dist(comb_index(i,:))).^2);
            %     position=fminsearch(cost_function,initial_guess);
            %     survive_posit(i, :) = position; 
            % end
            % x_mean = mean(survive_posit(:,1), 1);
            % y_mean = mean(survive_posit(:,2), 1);
            % z_mean = mean(survive_posit(:,3), 1);
            % x_std = std(survive_posit(:,1),1);
            % y_std = std(survive_posit(:,2),1);
            % z_std = std(survive_posit(:,3),1);
            % %position = [x_mean, y_mean, z_mean];
            % position(1) = mean(survive_posit((survive_posit(:,1)<=(x_mean+x_std*alpha)) & (survive_posit(:,1)>=(x_mean-x_std*alpha)), 1), 1);
            % position(2) = mean(survive_posit((survive_posit(:,2)<=(y_mean+y_std*alpha)) & (survive_posit(:,2)>=(y_mean-y_std*alpha)), 2), 1);
            % position(3) = mean(survive_posit((survive_posit(:,3)<=(z_mean+z_std*alpha)) & (survive_posit(:,3)>=(z_mean-z_std*alpha)), 3), 1);

        end

        %dist_fun=@(pos) sqrt(sum((bsxfun(@minus,XYZ_beacons,pos(:)')).^2,2));%
        %cost_function=@(pos) sum((dist_fun(pos)-measured_dist(:)).^2);%

        %position=fminsearch(cost_function,initial_guess);%
           

        record_posit = [position; record_posit(1:record_length-1,:)];
        record_posit
        Increase = 10;
        newX = linspace(1,20, 25 * Increase);
        record_posit_smooth = spline(x, y, newX);
        record_posit_smooth
        record_posit_smooth = record_posit_smooth(1:10,:);
        record_posit_smooth

        set(p9, 'XData', record_posit_smooth(:,1), 'YData', record_posit_smooth(:,2), 'ZData', record_posit_smooth(:,3))
        %set(p9, 'XData', position(1), 'YData', position(2), 'ZData', position(3))%
        formatSpec = '%.2f';
        th9.String = strcat("threshold: ", num2str(sum(index, "all")), " x : ",num2str(position(1),formatSpec)," y : ",num2str(position(2),formatSpec)," z : ",num2str(position(3),formatSpec));

        drawnow;
        pause(0.001);
    end
    
catch me  % if the error detected, close the radar
        
    radar1.close();
    
    clear radar1 frame1

    radar2.close();
    
    clear radar2 frame2

    radar3.close();
    
    clear radar3 frame3
        
    radar4.close();
    
    clear radar4 frame4

    radar5.close();
    
    clear radar5 frame5

    radar6.close();
    
    clear radar6 frame6

    radar7.close();
    
    clear radar7 frame7
        
    radar8.close();
    
    clear radar8 frame8
end

