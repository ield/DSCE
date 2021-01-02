fsa = 48000;    % Audio sampling rate
M = 5;          % Oversampling for radio bandwidth
fs = M*fsa;     % Radio sampling rate
fc= 1e9;      % Carrier frequency
Na = 4000;      % Length of audio buffer
Ny = Na*M;      % Length of receive buffer

m=1;            % AM modulation index (0 to 1)

chan = app.chan;                      % base-band channel
fc1=fsa*chan;
c = exp(1i*2*pi*(0:Ny-1)'*fc1/fs);    % Channel carrier

f0 = fsa/50;
xtest_a = sin(2*pi*(0:Na-1)'*f0/fsa);     % Test signal to be transmitted at audio rate
xtest_b = sin(2*pi*(0:Ny-1)'*f0/fs);      % Test signal to be transmitted at base-band rate

% Receive band-pass channel selection filter

Lr=128;
hr=fir1(Lr,(0.95*[-fsa, fsa]/2+fc1)*2/fs);
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
tx.Gain = -20;        % Gain goes from -89.75 dB to 0 dB
tx.BasebandSampleRate = fs;
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
disp('Radio on: Rx');

drawnow;
state = app.state;  % State of the radio. 0: OFF. 1: Rx. -1:Tx
state_0 = state;    % Preceding state 
chan_0 = chan;      % Preceding channel

while(state~=0)
%     x = ones(Ny, 1);
%     tx(x);  % Transmission buffer
    switch state
        case -1
            if state~=state_0; disp('Tx'); end
            
%             xa = adr();                                 % (1)
            xa = xtest_a;
            
            xa_upsample = upsample(xa, M);              % Upsampling
            
            [xb, zft] = filter(ht, 1, xa_upsample, zft);% Filter

            xc = xb*m/2 + (1 - m/2);                      % Scaling and sum, to be between 0 and 1.
%             xc = ones(Ny, 1);

            x = xc .* c;                                % Scale to channel a or b
%             x = ones(Ny, 1);
            tx(x);  % Transmission buffer
            
            
        otherwise % Commented so that tx and rx are working at the same
%         time when ptt is pressed.
            if state~=state_0; disp('Rx');  end
            
            % *** YOUR CODE ***
     
            y1=rx();                                    % Store buffers   
            y1 = double(y1)/2048;
            
            [y, zfr] = filter(hr, 1, y1, zfr);
            
            yc = abs(y);                                % Non-coherent demodulator
            
            yb = yc - mean(yc);                         % Delete the dc component
               
            ya = downsample(yb, M);                     % Downsample
            adw(ya);                                    % (1)
            
    end
    state_0 = state;
    drawnow;
    state = app.state;
    chan = app.chan;
    if chan~=chan_0
        disp(['Switched to channel ',num2str(chan)]);
        % *** YOUR CODE ***
        fc1=fsa*chan_0;
        c = exp(1i*2*pi*(0:Ny-1)'*fc1/fs);    % Channel carrier

        hr=fir1(Lr,(0.95*[-fsa, fsa]/2+fc1)*2/fs);
        zfr=zeros(Lr,1); % Filter memory

        ht=M*fir1(Lt,1/M);
        zft=zeros(Lt,1); % Filter memory
        chan_0 = chan;

    end
end

release(adw);
release(adr);

% Release transmitter and receiver
% *** YOUR CODE ***


disp('Radio off');

