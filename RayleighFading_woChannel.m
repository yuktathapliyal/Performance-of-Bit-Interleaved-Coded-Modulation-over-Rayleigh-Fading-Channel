% Communication system with QPSK mapping, Rayleigh fading channel, without channel coding
%==============================================================================
close all
clear all
clc
ModOrder = 4; % QPSK => Modulation Order = 4
BitRate = 10^6; % Transmission bit rate = Rb
fc = 10*10^9; % Carrier frequency
velocity = (60*10^3)/(3600); % v = 60km/h = (50/3) m/s
velocityOfLight = 3*10^8;
wavelength = velocityOfLight/fc;
fd = velocity/wavelength; % fd is the maximum Doppler shift
Tsig = 1/(BitRate/log2(ModOrder)); % Tsignal is duration of symbol after modulation and hence, after encoding as well.
Tcoh = 9/(16*pi*fd); % Tcoh is coherence time
NsymbolsAffected = ceil(Tcoh/Tsig); % Number of symbols affected by same fading coefficient (aplhai)
InputDataSize = 162*10^4; % Data bitstream length size
Nchannelcoeff = (InputDataSize/log2(ModOrder))/NsymbolsAffected; % Number of fading coefficients required for the modulated symbols stream.
EbNoVec = 0:5:30; % Eb/N0 in dB (not linear)
BER = zeros(1, round(length(EbNoVec))); % Pre-allocating BER vector for fast execution only.
for i = 1:length(EbNoVec)
  No = 10^(-1 * EbNoVec(i)/10);
  datastream = randi([0 1], [1 InputDataSize]);
  qpskmodsymbols = qpskmapping(datastream); % Gray-coded QPSK mapping
  alphavector = sqrt(0.5)*(randn(Nchannelcoeff) + 1i*randn(Nchannelcoeff)); % alphavector stores all the fading coefficients required.
  alphaMagnitude = abs(alphavector);
  receivedSymbols = zeros(1,round(length(qpskmodsymbols))); % Pre-allocation of vector for improving execution speed in the for loop.
  for indexalpha = 1:length(alphaMagnitude)
    FadingGroupedSymbols = qpskmodsymbols(1,((indexalpha-1)*NsymbolsAffected)+1:(indexalpha)*NsymbolsAffected);
    FadingPartialSignal = alphaMagnitude(1,indexalpha)*FadingGroupedSymbols; % This is the (alpha*x) in y=alpha*x + AWGnoise. But here it is part of whole x.
    AGWNoise = sqrt(No/2)*(randn(size(FadingGroupedSymbols)) + 1i * randn(size(FadingGroupedSymbols))); % Generating AWGN for each group of NsymbolsAffected with same alpha.
    receivedGroupedSymbols = FadingPartialSignal + AGWNoise; % Passing the encoded message through AWGN channel.
    EqualizedReceivedSymbols = receivedGroupedSymbols / alphaMagnitude(1,indexalpha);
    receivedSymbols(1,((indexalpha-1)*NsymbolsAffected+1:indexalpha*NsymbolsAffected)) = EqualizedReceivedSymbols; % Whole signal which is the combined version all the partial signals (grouped symbols).
  end
  qpskdemodsymbols = qpskdemapping(receivedSymbols); % Performinng demapping.
  BitErrorVec = (datastream ~= qpskdemodsymbols);
  BER(i) = sum(BitErrorVec)/length(BitErrorVec); % Creating a vector of BER of the signal
end
EbNoVecTheoretical = 0:30; % Eb/N0 % FOR THEORETICAL PLOT
BERuncodedTheoretical = 0.5 .* (1 - sqrt((10.^(EbNoVecTheoretical/10))./(1 + (10.^(EbNoVecTheoretical/10))))); % Theoretical uncoded BER for QPSK in fast and slow Rayleigh fading channel

%%PLOTS
figure (2)
semilogy(EbNoVec, BER, '-r+','linewidth',1.0)
axis([0.0 35 0 1.0])
xlabel('\it E_b/N_0 \rm(dB)', 'FontName', 'Times New Roman')
ylabel('Bit Error Rate', 'FontName', 'Times New Roman')
grid on
hold on
semilogy(EbNoVecTheoretical,BERuncodedTheoretical,'-b','linewidth',1.0)
legend('Scheme 1(a): Simulated without coding', 'Scheme 1(b): Theoretical without coding')
hold off

% FUNCTIONS
function qpskmodsymbols = qpskmapping(codeword) % QPSK Modulation function
  y = (2 * codeword) - 1;
  real = y(1:2:end);
  imag = y(2:2:end);
  qpskmodsymbols = real + (1i*imag);
end

function qpskdemodsymbols = qpskdemapping(receivedSymbols) % QPSK Demodulation function
  count = 1;
  modulOrder = 4; % Modulation order for QPSK = 4.
  demappedBitVec = zeros(1, length(receivedSymbols)*log2(modulOrder)); % Pre-allocating its size for fast execution of the program only. Other than this, it is not really required.
  for i=1:length(receivedSymbols)
    if real(receivedSymbols(i)) > 0
      demappedBitVec(count) = 1;
    else
      demappedBitVec(count) = 0;
    end
    if imag(receivedSymbols(i)) > 0
      demappedBitVec(count + 1) = 1;
    else
      demappedBitVec(count + 1) = 0;
    end
    count = count + 2;
  end
  qpskdemodsymbols=demappedBitVec;
end