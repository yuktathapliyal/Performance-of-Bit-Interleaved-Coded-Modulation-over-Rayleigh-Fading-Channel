% Communication system with QPSK mapping, Rayleigh fading channel,
(15, 11) Hamming Code (for channel coding)
%==============================================================================close all
clear all
clc
N = 15; % (15, 11) Hamming code, single error correcting code => t=1
K = 11;
M = N-K;
P = [1 1 1 1; % P = parity sub-matrix of generator matrix
    0 1 1 1; % P = [K x M]
    1 0 1 1;
    1 1 0 1;
    1 1 1 0;
    0 0 1 1;
    0 1 0 1;
    0 1 1 0;
    1 0 1 0;
    1 0 0 1;
    1 1 0 0];
G = [P eye(K)]; % generator matrix
H = [eye(M) P']; % parity check matrix, H = [I | P'], where P' is the transpose of P)
t = 1; % t = 1 as given
ModOrder = 4; % QPSK => Modulation Order = 4
BitRate = 10^6; % Transmission bit rate = Rb
fc = 10*10^9; % Carrier frequency
velocity = (60*10^3)/(3600); % v = 60km/h = (50/3) m/s
velocityOfLight = 3*10^8;
wavelength = velocityOfLight/fc;
fd = velocity/wavelength; % fd is the maximum Doppler shift
Tsig = 1/((BitRate*(N/K))/log2(ModOrder)); % Tsignal is duration of symbol after modulation and hence, after encoding as well.
Tcoh = 9/(16*pi*fd); % Tcoh is coherence time
NsymbolsAffected = ceil(Tcoh/Tsig); % Number of symbols affected by same fading coefficient (aplhai)
InputDataSize = 11*22*10^4; % Input data bitstream. It can be different from Rb.
EbNoVec = 0:5:35; % Eb/N0 in dB (not linear)
BER = zeros(1, round(length(EbNoVec))); % Pre-allocating BER vector for fast execution only.
datastream = randi([0 1], [1 InputDataSize]); % message stream
for i = 1:length(EbNoVec)
  No = 10^(-1 * EbNoVec(i)/10) * (N/K); % Computing Noise power with Eb normalized (fixed)
  codedstream = zeros(1, round(length(datastream)*(N/K))); % Pre-allocation of vector for improving execution speed in the for loop.
  flag = 0;
  startframe = 0; % index for framing
  for indexMsgFrame = 1:(length(datastream)/K) % index for moving frame-wise for message bits
    endframecount = startframe + K;
    message = datastream(startframe+1:endframecount); % framing stream by K=11 bits
    startframe = endframecount;
    codeword = mod((message * G), 2); % Encoding.
    codedstream(1,((indexMsgFrame-1)*N+1:indexMsgFrame*N)) = codeword; % Creating a vector of the whole coded stream of message bits.
  end
  if(mod(((length(codedstream)/log2(ModOrder))/NsymbolsAffected),2) ~= 0) % if symbols stream (signal) is not multiple of NsymbolsAffected (after encoding and modulation), then we pad zeros (dummy bits).
    flag = 1;
    Ndummysymbols = NsymbolsAffected - mod((length(codedstream)/log2(ModOrder)), NsymbolsAffected);
    Ndummybits = log2(ModOrder)* Ndummysymbols;
    codedstream = [codedstream zeros(1, Ndummybits)];
  end
  qpskmodsymbols = qpskmapping(codedstream); % Gray-coded QPSK mapping
  Nchannelcoeff = length(qpskmodsymbols)/NsymbolsAffected; % Number of channel (fading) coefficients that will be required to cater for all groups of 'NsymbolsAffected' symbols which will be affected by same alpha.
  alphavector = sqrt(0.5)*(randn(Nchannelcoeff) + 1i*randn(Nchannelcoeff)); % alphavector stores all the fading coefficients (alphai) required in complex form.
  alphaMagnitude = abs(alphavector);
  receivedSymbols = zeros(1,round(length(qpskmodsymbols))); % Pre-allocation of vector for improving execution speed in the for loop.
  for indexalpha = 1:length(alphaMagnitude) % index for moving alphai-wise.
    FadingGroupedSymbols = qpskmodsymbols(1,((indexalpha-1)*NsymbolsAffected)+1:(indexalpha)*NsymbolsAffected);
    FadingPartialSignal = alphaMagnitude(1,indexalpha)*FadingGroupedSymbols; % This is the (alpha*x) in y = alpha*x + AWGnoise. But here it is part of the whole x.
    AGWNoise = sqrt(No/2)*(randn(size(FadingGroupedSymbols)) + 1i * randn(size(FadingGroupedSymbols))); % Generating AWGN for each group of NsymbolsAffected with same alpha.
    receivedGroupedSymbols = FadingPartialSignal + AGWNoise; % Passing the encoded message through AWGN channel.
    EqualizedReceivedSymbols = receivedGroupedSymbols / alphaMagnitude(1,indexalpha);
    receivedSymbols(1,((indexalpha-1)*NsymbolsAffected+1:indexalpha*NsymbolsAffected)) = EqualizedReceivedSymbols; % Creating a vector for the whole signal which is the combined version all the partial signals (grouped symbols).
  end
  qpskdemodsymbols = qpskdemapping(receivedSymbols); % Performinng demapping. qpskdemodsymbols includes the dummy bits, if they were added before. Hence need to remove them, if any.
  if(flag == 1) % If flag==1, then it means dummy bits were added before.
    qpskdemodsymbols(end+1-Ndummybits:end) = []; % Discarding the dummy bits
  end
  Ie = eye(N); % Rows of Ie are error patterns (vectors)
  syndromeTable = mod((Ie * double(H')), 2); % Corresponding syndrome for each single-bit error pattern/vector
  startframe = 0; % index for start of frame
  clear endframecount;
  BERperFrameVec = zeros(1,(length(qpskdemodsymbols)/N)); % Pre-allocating its size with zeros for fast execution of the program only.
  for indexRXcodwrd = 1:(length(qpskdemodsymbols)/N) % index for moving codeword-wise in the long demapped stream
    endframecount = startframe + N;
    demappedCodeword = qpskdemodsymbols(startframe+1:endframecount); % framing stream by K=11 bits
    syndrome = mod((demappedCodeword * double(H')), 2); % Syndrome
    if (all(syndrome == 0))
      decodedvec = demappedCodeword; % if S = 0, then demodulated vectors is accepted as codeword
    else
      lookSyndromeinIe = ismember(syndromeTable, syndrome, 'rows'); % else look up in syndrome table for match
      if (all(lookSyndromeinIe == 0)) % All-zero means that no match is found.
        decodedvec = demappedCodeword;
      else
        pidx = find(lookSyndromeinIe == 1); % If a match is found
        decodedvec = mod((demappedCodeword + Ie(pidx, :)),2); % Adding the corresponding e to demapped vector
      end
    end
    decodedstream((indexRXcodwrd-1)*K+1: indexRXcodwrd*K)=decodedvec(M+1:N);
    startframe = endframecount;
  end
  BER(i) = sum(mod((decodedstream + datastream),2))/length(datastream);
end

EbNoVecTheoretical = 0:35; % Eb/N0 % FOR UN-CODED THEORETICAL PLOT
BERuncodedTheoretical = 0.5 .* (1 - sqrt((10.^(EbNoVecTheoretical/10))./(1 + (10.^(EbNoVecTheoretical/10))))); % Theoretical uncoded BER for QPSK in fast and slow Rayleigh fading channel
%%PLOTS
figure (2)
semilogy(EbNoVec, BER, '-r+','linewidth',1.0)
axis([0.0 55 10^-8 1.0])
xlabel('\it E_b/N_0 \rm(dB)', 'FontName', 'Times New Roman')
ylabel('Bit Error Rate', 'FontName', 'Times New Roman')
grid on
hold on
semilogy(EbNoVecTheoretical,BERuncodedTheoretical,'-b','linewidth',1.0)
legend('Scheme 2(a): Simulated with Hamming code', 'Scheme 1(b): Theoretical without coding')
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
  demappedBitVec = zeros(1, length(receivedSymbols)*log2(modulOrder)); % Pre-allocating its size for fast execution of the program only.
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
    
