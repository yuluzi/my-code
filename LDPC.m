mod = comm.QPSKModulator('BitInput',true); % ����������
demod = comm.QPSKDemodulator('BitOutput',true,...
  'DecisionMethod','Approximate log-likelihood ratio', ...
  'VarianceSource','Input port'); % ���������
chan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)'); % AWGN�ŵ�
err = comm.ErrorRate;
enc = comm.LDPCEncoder; % ����LDPC������
dec = comm.LDPCDecoder; % ����LDPC������
snrVec = [0 0.2 0.4 0.6 0.65 0.7 0.75 0.8]; % ����ȷ�Χ
ber = zeros(length(snrVec),1);
for k = 1:length(snrVec)
  chan.SNR = snrVec(k);
  noiseVar = 1/10.^(snrVec(k)/10);
  errorStats = zeros(1,3);
  while errorStats(2) <= 200 && errorStats(3) < 5e6
      data = logical(randi([0 1],32400,1));     % ��������������
      encData = step(enc,data);              % LDPC����
      modSig = step(mod,encData);          % ����
      rxSig = step(chan,modSig);            % AWGN�ŵ�
      demodSig = step(demod,rxSig,noiseVar); % ���
      rxData = step(dec,demodSig);          % LDPC����
      errorStats = step(err,data,rxData);       % ͳ�ƴ���
  end
  ber(k) = errorStats(1);
  reset(err)
end
semilogy(snrVec,ber,'-o')
grid
xlabel('SNR (dB)')
ylabel('�������')