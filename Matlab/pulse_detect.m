% Setup
clear; clc; close all;

% Create pulse to detect

rng('default');%random variable will be selected using default algorithm
PulseLen = 64;%64 bits will be the length of the pulse
theta = rand(PulseLen,1);%an array is generated with size 64x1 of values 0 to 1
pulse = exp(1i*2*pi*theta);%the values are converted into complex form

% Insert pulse to Tx signal

rng('shuffle');%random variable will be selected using default algorithm
TxLen = 5000;%total transmition length will be 5000 bits
PulseLoc = randi(TxLen-PulseLen*2);% defining where to insert the pulse

TxSignal = zeros(TxLen,1);% the rest 5000-64 = 4936 bits should be zero while transmitting
TxSignal(PulseLoc:PulseLoc+PulseLen-1) = pulse; %merging both pulse to 5000 bits chain

% Create Rx signal by adding noise

Noise = complex(randn(TxLen,1),randn(TxLen,1));%random complex numbers
RxSignal = TxSignal + Noise; % addding noise to signal

% Scale Rx signal to +/- one

scale1 = max([abs(real(RxSignal)); abs(imag(RxSignal))]); %converting it for convinience
RxSignal = RxSignal/scale1;%dividing it by max value of signal

% Create matched filter coefficients

CorrFilter = conj(flip(pulse))/PulseLen;%step 1 of corelation filter

% Correlate Rx signal against matched filter

FilterOut = filter(CorrFilter,1,RxSignal);%step 2 of corelation filter

% Find peak magnitude & location

[peak, location] = max(abs(FilterOut));%getting the value of peak at filter output

% Print results
figure(1)%creating a figure
subplot(311); plot(real(TxSignal)); title('Tx Signal (real)');%the signal without noise
subplot(312); plot(real(RxSignal)); title('Rx Signal (real)');%the signal with noise

t = 1:length(FilterOut);%length of signal
str = sprintf('Peak found at %d with a value of %.3d',location,peak);%print the value of peak
subplot(313); plot(t,abs(FilterOut),location,peak,'o'); title(str);%the filtered signal

% Implementing the same above mentioned algorithm in a hardware friendly manner.

WindowLen = 11;%comparing the 5000 bits with 11 bits at a time serially
MidIdx = ceil(WindowLen/2);%gets the value of 6th bit in the series
threshold = 0.03;%only the peak would be more than this value.

%note: for convenience the threshold value was adjusted from experimental
       %trial and error method

% Compute magnitude squared to avoid sqrt operation
MagSqOut = abs(FilterOut).^2;

% Sliding window operation
for n = 1:length(FilterOut)-WindowLen %1 to 4989

    % Compare each value in the window to the middle sample via subtraction
    DataBuff = MagSqOut(n:n+WindowLen-1);%stream of 11 bits
    MidSample = DataBuff(MidIdx);%bit 6 of the sample
    CompareOut = DataBuff - MidSample; % this is a vector

    % if all values in the result are negative and the middle sample is
    % greater than a threshold, it is a local max
    if all(CompareOut <= 0) && (MidSample > threshold)
        peak_2 = MidSample;
        location_2 = n + (MidIdx-1);%gives the peak location
    end
end

% Simulate model in fixed-point or floating-point
fxpt_mode = true;
if fxpt_mode  % fixed-point
    DT_input = fixdt(1,16,14);%fixdt(signed/unsigned',wordlength,fractionlength)
    DT_filter = fixdt(1,18,15);
    DT_power = fixdt(1,18,11);
else  % floating-point
    DT_input = 'double';%64 bit floating point
    DT_filter = 'double';
    DT_power = 'double';
end
DT_coeff = fixdt(1,18); % coeff is treated as double if input is double

if iscolumn(CorrFilter)
    CorrFilter = transpose(CorrFilter); % need row vector for filter block
end

%stimulation time to be length of the signal plus 30 time units
%30 time units= 10 in front end + 10 in back end + 10 buffer
SimTime = length(RxSignal) + WindowLen + 30;

% Simulate model
slout = sim('project_detect');

% Correlation filter output
FilterOutSL = getLogged(slout,'filter_out');
FilterValid = getLogged(slout,'filter_valid');
FilterOutSL = FilterOutSL(FilterValid);%signal through filter

%to comapare and print logged signal values
compareData(real(FilterOut),real(FilterOutSL),{2 3 1},'ML vs SL correlator output (re)');
compareData(imag(FilterOut),imag(FilterOutSL),{2 3 2},'ML vs SL correlator output (im)');

% Magnitude squared output
MagSqSL = getLogged(slout,'mag_sq_out');
MagSqSL = MagSqSL(FilterValid);%magnitude square of filter

%to comapare and print logged signal values
compareData(MagSqOut,MagSqSL,{2 3 3},'ML vs SL mag-squared output');

% Peak value
MidSampleSL = getLogged(slout,'mid_sample');
Detected = getLogged(slout,'detected');
PeakSL = double(MidSampleSL(Detected>0));

fprintf('\nPeak location = %d, magnitude = %.3d using global max\n',location,peak);
fprintf('Peak location = %d, mag-squared = %.3d using local max\n',location_2,peak_2);
fprintf('Peak mag-squared from Simulink = %.3d, error = %.3d\n',PeakSL,abs(peak_2-PeakSL));

%%Defining the functions used in the code above

%getLogged function
function signal_val = getLogged(simout_obj,signal_name)%function instantiation

logsout = simout_obj.logsout;%load logged signal 

%checking for possible error in outputs
if isempty(logsout)
    error('No logged signal found. Make sure ''%s'' is logged in the model',...
        signal_name);
end

sig = logsout.getElement(signal_name);%read elements of logged signal

%checking for possible error in outputs
if isempty(sig)
    error('Signal ''%s'' not found. Make sure it is logged and named correctly.',...
        signal_name);
end

signal_val = squeeze(sig.Values.Data);
%squeeze returns an array with the same elements as the input,...
%but with dimensions of length 1 removed.

end

%compareData function
%function instantiation
function err_vec = compareData(reference,actual,figure_number,textstring)

% Vector input only
if ~isvector(reference) || ~isvector(actual)
    error('Input signals must be vector');%check for inputs to be vector
else
    if isrow(reference)
        reference = transpose(reference);%convert to column vector
    end
    if isrow(actual)
        actual = transpose(actual);%convert to column vector
    end
end

% Make signals same length if necessary
if length(reference) ~= length(actual)
%     warning(['Length of reference (%d) is not the same as actual signal (%d).'...
%         ' Truncating the longer input.'],length(reference),length(actual));
    
    %balance length of both the vectors to be equal
    len = 1:min(length(reference),length(actual));
    reference = reference(len);%length modification of reference
    actual = actual(len);%length modification of actual
end

% Turn complex into vector
if xor(isreal(reference),isreal(actual))%checking for nature of both signals
    error('Input signals are not both real or both complex');
elseif ~isreal(reference)
    ref_vec = double([real(reference) imag(reference)]);%covert to double
    act_vec = double([real(actual) imag(actual)]);%covert to double
    tag = {'(Real)','(Imag)'};
else
    ref_vec = double(reference);%covert to double
    act_vec = double(actual);%covert to double
    tag = {''};
end

% Configure figure
if iscell(figure_number)
    if size(ref_vec,2) > 1 % complex
        error('Cannot yet subplot multiple complex inputs');
    else
        figure(figure_number{1})%defines figure(2)
    end
else
    figure(figure_number)
end
c = get(groot,'defaultAxesColorOrder');%used at plotting

% Compute error
err_vec = ref_vec - act_vec;%error vector to be printed on output
max_err = max(abs(err_vec));%max error to be printed on output
max_ref = max(abs(ref_vec));%max reference to be printed on output

%output printed on command window
fprintf('\nMaximum error for %s out of %d values\n',textstring,length(actual));

for n = 1:size(ref_vec,2)
   
    %output printed on Command window
    fprintf('%s %d (absolute), %d (percentage)\n',tag{n},max_err(n),max_err(n)/max_ref(n)*100);
   
    if iscell(figure_number)
        total_plot = figure_number{2};% equal to 3
        plot_num = figure_number{3};%figure number can be 1 to 3
   
    else
        total_plot = size(ref_vec,2);
        plot_num = n;
    end
    
    subplot(total_plot,1,plot_num)%subplot generation
    
    plot(ref_vec(:,n),'Color',c(3,:));%uses 1st color for reference
    
    hold on %plots 3 graphs over the same plot
    
    plot(act_vec(:,n),'Color',c(1,:));%uses 2nd color for actual
    plot(err_vec(:,n),'Color',c(2,:));%uses 3rd color for error
    
    legend('Reference','Actual','Error')%defining legend
    hold off 
    
    title(sprintf('%s %s, max error = %.3d',textstring,tag{n},max_err(n)));%title
end
end
