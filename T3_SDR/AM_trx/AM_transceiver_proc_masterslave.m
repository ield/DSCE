fsa = 48000;    % Audio sampling rate
M = 5;          % Oversampling for radio bandwidth
fs = M*fsa;     % Radio sampling rate
fc= 1e9;      % Carrier frequency
Na = 4000;      % Length of audio buffer
Ny = Na*M;      % Length of receive buffer

m=1;            % AM modulation index (0 to 1)

% The initial configuration is the slave configuration

chan_tx = 2;                      % base-band channel
fc1_tx=fsa*chan_tx;
c_tx = exp(1i*2*pi*(0:Ny-1)'*fc1_tx/fs);    % Channel carrier

chan_rx = 1;
fc1_rx=fsa*chan_rx;

f0 = fsa/50;
xtest_a = sin(2*pi*(0:Na-1)'*f0/fsa);     % Test signal to be transmitted at audio rate
xtest_b = sin(2*pi*(0:Ny-1)'*f0/fs);      % Test signal to be transmitted at base-band rate

% Receive band-pass channel selection filter

Lr=128;
hr=fir1(Lr,(0.95*[-fsa, fsa]/2+fc1_rx)*2/fs);
zfr=zeros(Lr,1); % Filter memory

% Transmit interpolation low-pass filter

Lt=128;
ht=M*fir1(Lt,1/M);
zft=zeros(Lt,1); % Filter memory

% Open audio reader
adr=audioDeviceReader(fsa,'SamplesPerFrame',Na);

% Open audio writer
adw=audioDeviceWriter(fsa,'BufferSize',Na);

%% Open transmitter and configure it. Release it at the end
tx=sdrtx('Pluto');

% Set transmitter parameters:
% Carrier frequency
% Sampling rate
% Transmitter gain. -10 is adequate for connecting the spectrum analyser

tx.CenterFrequency = fc;
tx.BasebandSampleRate = fs;
tx.Gain = 0;        % Gain goes from -89.75 dB to 0 dB
% tx.FrequencyCorrection = 5.45;     %Nacho
% tx.FrequencyCorrection = 2;     %Diego

%% Open receiver and configure it. Release it at the end
rx=sdrrx('Pluto');

% Set receiver parameters:
% Carrier frequency
% Sampling rate
% AGC
% Receiver gain in dB

rx.CenterFrequency = fc;
rx.BasebandSampleRate = fs;
rx.GainSource = 'Manual';       % Disable AGC for manual setting of gain
rx.Gain = 50;                    % The gain goes from -4 to 62 dB

%%
disp('Radio on: Slave');

drawnow;
state = app.state;  % State of the radio. 0: OFF. 1: Rx. -1:Tx
state_0 = state;    % Preceding state 
% chan_0 = chan;      % Preceding channel

while(state~=0)        
    xa = adr();                                 % (1)

    xa_upsample = upsample(xa, M);              % Upsampling

    [xb, zft] = filter(ht, 1, xa_upsample, zft);% Filter
%             xb = xa_upsample;
    xc = xb*m/2 + 1 - m/2;                      % Scaling and sum, to be between 0 and 1.

    x = xc .* c_tx;                                % Scale to channel a or b

    tx(x);  % Transmission buffer

    y1=rx();                                    % Store buffers   
    y1 = double(y1)/2048;

    [y, zfr] = filter(hr, 1, y1, zfr);

    yc = abs(y);                                % Non-coherent demodulator

    yb = yc - mean(yc);                         % Delete the dc component

    ya = downsample(yb, M);                     % Downsample
    adw(ya);                                    % (1)
            
    state = app.state;
    drawnow;
    
    % When it changes between slave and master: it is like pressing the 
    % ptt, so initially, all radios are slave. What it is done is to
    % change the carriers and filter of the slave and master: a change in
    % variables.
    if state~=state_0      
        
        chan_tx_old = chan_tx;                      % base-band channel
        chan_tx = chan_rx;
        chan_rx = chan_tx_old;
        
        fc1_tx=fsa*chan_tx;
        c_tx = exp(1i*2*pi*(0:Ny-1)'*fc1_tx/fs);    % Channel carrier

        fc1_rx=fsa*chan_rx;
        hr=fir1(Lr,(0.95*[-fsa, fsa]/2+fc1_rx)*2/fs);
        zfr=zeros(Lr,1); % Filter memory
        freqz(hr, 1, 'whole');

        % Transmit interpolation low-pass filter
        ht=M*fir1(Lt,1/M);
        zft=zeros(Lt,1); % Filter memory
        
        %5. Refresh the state and inform of the new state
        state_0 = state;
        if state_0 == -1
            disp('The radio is master');
        elseif state_0 == 1
            disp('The radio is slave');
        end
    end
end

release(adw);
release(adr);

% Release transmitter and receiver
% *** YOUR CODE ***


disp('Radio off');

