mod = comm.QPSKModulator('BitInput',true); % 创建调制器
demod = comm.QPSKDemodulator('BitOutput',true,...
  'DecisionMethod','Approximate log-likelihood ratio', ...
  'VarianceSource','Input port'); % 创建解调器
chan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)'); % AWGN信道
err = comm.ErrorRate;
enc = comm.LDPCEncoder; % 创建LDPC编码器
dec = comm.LDPCDecoder; % 创建LDPC译码器
snrVec = [0 0.2 0.4 0.6 0.65 0.7 0.75 0.8]; % 信噪比范围
ber = zeros(length(snrVec),1);
for k = 1:length(snrVec)
  chan.SNR = snrVec(k);
  noiseVar = 1/10.^(snrVec(k)/10);
  errorStats = zeros(1,3);
  while errorStats(2) <= 200 && errorStats(3) < 5e6
      data = logical(randi([0 1],32400,1));     % 产生二进制数据
      encData = step(enc,data);              % LDPC编码
      modSig = step(mod,encData);          % 调制
      rxSig = step(chan,modSig);            % AWGN信道
      demodSig = step(demod,rxSig,noiseVar); % 解调
      rxData = step(dec,demodSig);          % LDPC译码
      errorStats = step(err,data,rxData);       % 统计错误
  end
  ber(k) = errorStats(1);
  reset(err)
end
semilogy(snrVec,ber,'-o')
grid
xlabel('SNR (dB)')
ylabel('误比特率')