M=16; % 调制阶数 调制方式16-QAM 常见的BPSK QPSK 16-QAM 64-QAM 256-QAM
EbNo=-5:-1; % 信噪比范围 Eb代表平均到每个比特上的信号能量，N0代表噪声的功率谱密度。
frmLen=500; % 帧长·                  
ber=zeros(size(EbNo));%初始化误码率矩阵
enc=comm.TurboEncoder('InterleaverIndicesSource','Input port'); % 创建编码器 comm是通信工具箱
dec=comm.TurboDecoder('InterleaverIndicesSource','Input port','NumIterations',3); % 创建译码器 迭代次数为4
mod=comm.RectangularQAMModulator('ModulationOrder',M, ...
  'BitInput',true,'NormalizationMethod','Average power'); % 创建QAM调制器
%我们要传输的信号一般都为低频弱电信号。要远距离发射必须要用天线，
%而天线的长度必须大于等于4倍的波长才能有效不失真传输。
%频率=光速/波长，因为低频的波天线就要做的非常长，那是不可能实现的。
%所以就要用低频的信号（你想要发射的)加载在高频载波上再传输，这个过程就叫做调制!
demod=comm.RectangularQAMDemodulator('ModulationOrder',M, ...% 创建QAM译制器
  'BitOutput',true,'NormalizationMethod','Average power', ...%归一化方法-平均功率
  'DecisionMethod','Log-likelihood ratio', ...%调制决策方法-LLR 默认Hard decision| Log-likelihood ratio| Approximate log-likelihood ratio. 
  'VarianceSource','Input port'); %噪声方差的来源
chan=comm.AWGNChannel('EbNo',EbNo,'BitsPerSymbol',log2(M)); %每个符号的位数 4
% 创建AWGN信道
errorRate = comm.ErrorRate;%创建一个计算错误率的目标，通过将接收数据与发射数据比较的方式得到错误率。
% 1.选择随机数据包长度，生成随机二进制数据。
% 2.计算输出码字长度和 编码率。
% 3.计算信噪比 (SNR) 和噪声方差。
% 4.生成交织器索引。
% 5.Turbo 编码数据。
% 6.应用 16-QAM 调制，并对平均信号功率进行归一化。
% 7.将调制信号通过 AWGN 通道。
% 8.使用LLR方法解调带噪信号，输出软比特，并对平均信号功率进行归一化。
% 9.Turbo 解码数据。由于来自解调器的位映射顺序与 Turbo 解码器预期的映射顺序相反。
% 10.计算误差统计量。
for k = 1:length(EbNo)
  errorStats = zeros(1,3); 
  noiseVar = 10^(-EbNo(k)/10)*(1/log2(M));
  chan.EbNo = EbNo(k);
  while errorStats(2) < 100 && errorStats(3) < 1e7
      data = randi([0 1],frmLen,1); % 创建二进制随机数据
      intrlvrInd = randperm(frmLen); % 交织
      %调用step方法，执行对象
      encodedData = step(enc,data,intrlvrInd); % Turbo编码
      modSignal = step(mod,encodedData); % 调制
      receivedSignal = step(chan,modSignal); % AWGN信道
      demodSignal = step(demod,receivedSignal,noiseVar); % 解调
      receivedBits = step(dec,-demodSignal,intrlvrInd); % 译码
      errorStats = step(errorRate,data,receivedBits); % 统计差错
  end
  %保存BER数据并重置误码率对象
  ber(k) = errorStats(1);
  reset(errorRate)
end
semilogy(EbNo,ber,'r-o')
grid %切换改变主网格线的可见性
xlabel('Eb/No (dB)')
ylabel('误比特率')
uncodedBER = berawgn(EbNo,'qam',M); % 未编码的BER
hold on %保持原图
semilogy(EbNo,uncodedBER,'b-')
legend('Turbo编码','未编码')%图例