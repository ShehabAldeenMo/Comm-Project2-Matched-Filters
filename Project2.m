%% Clear Section---------------------------------
% Clear the workspace, command window, and close all figures.
% This section ensures a clean slate before running the simulations.
clear; clc; close all;

%% Random Signal Generation -------------------------------------------
% Generate a random sequence of 10 bits
randomBits = randi([0, 1], 1, 10);
randomBitsLarge = randi([0, 1], 1, 10000);

% Convert bits to bipolar (-1, +1)
bipolarSequence = (2 * randomBits - 1);
bipolarSequenceLarge = (2 * randomBitsLarge - 1);

% Upsample the sequence to sample every 200 ms
sampledSignal = upsample(bipolarSequence, 5);
sampledSignalLarge = upsample(bipolarSequenceLarge, 5); 

% Define the pulse shaping function
pulseShape = [5, 4, 3, 2, 1] / sqrt(55);

% Transmit the signal through the pulse shaping filter
transmittedSignal = conv(sampledSignal, pulseShape);
transmittedSignalLarge = conv(sampledSignalLarge, pulseShape);

% Remove zero padding introduced by convolution
transmittedSignal = transmittedSignal(:, 1:50);
transmittedSignalLarge = transmittedSignalLarge(:, 1:50000);

%% Matched Filter -----------------------------------------------------
% Generate the matched filter by flipping the pulse shaping function
matchedFilter = fliplr(pulseShape);

% Apply the matched filter to the transmitted signal
matchedFilterOutput = conv(matchedFilter, transmittedSignal);

% Sample the matched filter output every second
matchedFilterOutputSampled = matchedFilterOutput(5:5:end);

%% Normalized Filter --------------------------------------------------
% Define the coefficients for the normalized filter
normalizedFilterCoefficients = [1, 1, 1, 1, 1] / sqrt(5);

% Apply the normalized filter to the transmitted signal
normalizedFilterOutput = conv(transmittedSignal, normalizedFilterCoefficients);

% Sample the normalized filter output every second
normalizedFilterOutputSampled = normalizedFilterOutput(5:5:end);

%% Plotting -----------------------------------------------------------
% This section creates various plots to visualize the filter outputs.

figure;
subplot(2, 1, 1);
plot(matchedFilterOutput, 'b', 'linewidth', 1);
xlabel("Time");
title("Output of the matched filter");
subplot(2, 1, 2);
plot(normalizedFilterOutput, 'r', 'linewidth', 1);
xlabel("Time");
title("Output of the normalized filter");

figure;
subplot(2, 1, 1);
stem(matchedFilterOutput, 'b', 'linewidth', 1);
xlabel("Time");
title("Output of the matched filter before sampling");
subplot(2, 1, 2);
stem(normalizedFilterOutput, 'r', 'linewidth', 1);
xlabel("Time");
title("Output of the normalized filter before sampling");

figure;
subplot(2, 1, 1);
stem(matchedFilterOutputSampled, 'b', 'linewidth', 1);
xlabel("Time");
title("Output of the matched filter after sampling");
subplot(2, 1, 2);
stem(normalizedFilterOutputSampled, 'r', 'linewidth', 1);
xlabel("Time");
title("Output of the normalized filter after sampling");

%% Correlator in a Noise-Free System ---------------------------------
% Initialize variables
correlatorInput = repmat(pulseShape, 1, 10) .* transmittedSignal;
windowSize = length(pulseShape);  % Size of the pulse shape

% Preallocate correlatorOutput to its maximum possible size
correlatorOutput = zeros(1, 50);

% Calculate correlator output
for i = 1:50
    startIndex = max(1, i - windowSize + 1);
    endIndex = i;
    correlatorOutput(i) = sum(correlatorInput(startIndex:endIndex));
end

% Sample the correlator output
correlatorOutputSampled = correlatorOutput(5:5:end);

%% Plotting -----------------------------------------------------------
% This section plots the correlator and matched filter outputs.
figure;
plot(correlatorOutput, 'r', 'linewidth', 1);
hold on;
plot(matchedFilterOutput, 'b', 'linewidth', 1);
legend('Correlator Output', 'Matched Filter Output');
title('Output of Correlator & Matched Filter');

% Plot correlator output and matched filter output without sampling
figure;
stem(correlatorOutput, 'r', 'linewidth', 1);
hold on;
stem(matchedFilterOutput, 'b', 'linewidth', 1);
legend('Correlator Output', 'Matched Filter Output');
title('Output of Correlator & Matched Filter before Sampling');

% Plot sampled correlator output and matched filter output
figure;
stem(correlatorOutputSampled, 'r', 'linewidth', 1);
hold on;
stem(matchedFilterOutputSampled, 'b', 'linewidth', 1);
legend('Correlator Output', 'Matched Filter Output');
title('Output of Correlator & Matched Filter after Sampling');


%% Noisy Environment---------------------------------
% Generate Gaussian noise
noise = randn(1, size(transmittedSignalLarge, 2)); 

% Define noise power and signal power
EbNo = 1;
Eb = 1;

% Pre-allocate arrays to store BER values
numIterations = length(-2:5);
EbNoRatio = zeros(1, numIterations);
MF_BER_values = zeros(1, numIterations);
NF_BER_values = zeros(1, numIterations);

% Iterate through different Eb/No ratios
for j = -2:5
    EbNo = Eb / (10^(j/10)); 
    noiseScaled = noise .* sqrt(EbNo / 2);
    
    % Add noise to signal
    noisySignal = transmittedSignalLarge + noiseScaled; 
    
    % Matched filter
    matchedFilter = fliplr(pulseShape); 
    matchedFilterOutput = conv(noisySignal, matchedFilter); 
    
    % Normalized filter
    normalizedFilter = [5 5 5 5 5] / sqrt(125); 
    normalizedFilterOutput = conv(noisySignal, normalizedFilter); 
    
    % Downsample the filtered signals
    downsampledMF = matchedFilterOutput(5:5:end); 
    downsampledNF = normalizedFilterOutput(5:5:end); 
    
    % Thresholding for matched filter
    downsampledMF(downsampledMF < 0) = -1; 
    downsampledMF(downsampledMF > 0) = 1; 
    
    % Thresholding for normalized filter
    downsampledNF(downsampledNF < 0) = -1; 
    downsampledNF(downsampledNF > 0) = 1; 
    
    % Calculate Bit Error Rate (BER)
    [numErrorsMF, MF_BER] = symerr(downsampledMF, bipolarSequenceLarge); 
    [numErrorsNF, NF_BER] = symerr(downsampledNF, bipolarSequenceLarge); 
    
    % Store BER values
    EbNoRatio(j+3) = Eb/EbNo; 
    MF_BER_values(j+3) = MF_BER; 
    NF_BER_values(j+3) = NF_BER; 
end

% Calculate theoretical BER
theoreticalBER = 0.5 * erfc(sqrt(EbNoRatio)); 

%% Plotting -----------------------------------------------------------
% This section plots the BER results.

figure; 
semilogy(-2:5, MF_BER_values, 'b', 'linewidth', 1); 
hold on; 
semilogy(-2:5, NF_BER_values, 'r', 'linewidth', 1); 
hold on; 
semilogy(-2:5, theoreticalBER, 'g', 'linewidth', 1); 
xlabel('Eb/No'); 
ylabel('BER'); 
legend('MF BER', 'NF BER', 'Theoretical BER'); 
title('Matched Filter BER, Normalized Filter BER & Theoretical BER');

%% ISI & Raised Cosine----------------------------------------------------
% Define roll-off factors [0 0 1 1]
rollOffFactors = [0, 0, 1, 1];

% Define delays [2 8 2 8]
delays = [2, 8, 2, 8];

% Generate random bit sequence
randomBits = randi([0, 1], 1, 100);

% Map the data from [0, 1] to [-1, 1]
mappedData = ((2 * randomBits) - 1) * 1;

% Introduce more data (zeros) to upsample the data to higher sampling rate
upsampledData = upsample(mappedData, 5);

% Loop through different roll-off factors and delays
for i = 1:4

    % Get the numerator and denominator of the raised cosine transfer function
    h = rcosdesign(rollOffFactors(i), delays(i) * 5, 5, 'sqrt');

    % Apply the transfer function to the data
    filteredTransmittedSignal = filter(h, 1, upsampledData);

    % Plot the eye pattern at transmitter
    eyediagram(filteredTransmittedSignal, 10);
    title(sprintf('Eye Diagram at Transmitter (Roll-off: %d, Delay: %d)', rollOffFactors(i), delays(i)));
    xlabel('Time');
    ylabel('Amplitude');
    figure;

    % Apply matched filter
    filteredReceivedSignal = filter(h, 1, filteredTransmittedSignal);

    % Plot the eye pattern at receiver
    eyediagram(filteredReceivedSignal, 10);
    title(sprintf('Eye Diagram at Receiver (Roll-off: %d, Delay: %d)', rollOffFactors(i), delays(i)));
    xlabel('Time');
    ylabel('Amplitude');
    figure;
end 