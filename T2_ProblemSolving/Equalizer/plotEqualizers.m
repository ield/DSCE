function [] = plotEqualizers(b, f, y, dec)
%Plots all the required data for equalizers
figure;
subplot(2, 3, 1);
stem(0:length(b)-1, b); % Channel in time domain

subplot(2, 3, 2);
stem(0:length(f)-1, f); % Equalizer in time domain

subplot(2, 3, 3);
prod = conv(b, f);
stem(0:length(prod)-1, prod);   %Result channel and equalizer time domain
% fprintf('e. The result in t must be the highest possible delta centered in the delay\n');

subplot(2, 3, 4);
plot(abs(fft(b))); % Channel in f domain

subplot(2, 3, 5);
plot(abs(fft(f))); % Equalizer in f domain

subplot(2, 3, 6);
plot(abs(fft(prod))); %Result channel and equalizer freq domain
% fprintf('The result in f must be close to a constant (the transformation of the delta)\n');

% Plot the eye diagram
figure;
subplot(1, 2, 1);
plot(y, 'o');

subplot(1, 2, 2);
plot(dec, 'o');
end

