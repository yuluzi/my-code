M=16; % ���ƽ��� ���Ʒ�ʽ16-QAM ������BPSK QPSK 16-QAM 64-QAM 256-QAM
EbNo=-5:-1; % ����ȷ�Χ Eb����ƽ����ÿ�������ϵ��ź�������N0���������Ĺ������ܶȡ�
frmLen=500; % ֡����                  
ber=zeros(size(EbNo));%��ʼ�������ʾ���
enc=comm.TurboEncoder('InterleaverIndicesSource','Input port'); % ���������� comm��ͨ�Ź�����
dec=comm.TurboDecoder('InterleaverIndicesSource','Input port','NumIterations',3); % ���������� ��������Ϊ4
mod=comm.RectangularQAMModulator('ModulationOrder',M, ...
  'BitInput',true,'NormalizationMethod','Average power'); % ����QAM������
%����Ҫ������ź�һ�㶼Ϊ��Ƶ�����źš�ҪԶ���뷢�����Ҫ�����ߣ�
%�����ߵĳ��ȱ�����ڵ���4���Ĳ���������Ч��ʧ�洫�䡣
%Ƶ��=����/��������Ϊ��Ƶ�Ĳ����߾�Ҫ���ķǳ��������ǲ�����ʵ�ֵġ�
%���Ծ�Ҫ�õ�Ƶ���źţ�����Ҫ�����)�����ڸ�Ƶ�ز����ٴ��䣬������̾ͽ�������!
demod=comm.RectangularQAMDemodulator('ModulationOrder',M, ...% ����QAM������
  'BitOutput',true,'NormalizationMethod','Average power', ...%��һ������-ƽ������
  'DecisionMethod','Log-likelihood ratio', ...%���ƾ��߷���-LLR Ĭ��Hard decision| Log-likelihood ratio| Approximate log-likelihood ratio. 
  'VarianceSource','Input port'); %�����������Դ
chan=comm.AWGNChannel('EbNo',EbNo,'BitsPerSymbol',log2(M)); %ÿ�����ŵ�λ�� 4
% ����AWGN�ŵ�
errorRate = comm.ErrorRate;%����һ����������ʵ�Ŀ�꣬ͨ�������������뷢�����ݱȽϵķ�ʽ�õ������ʡ�
% 1.ѡ��������ݰ����ȣ�����������������ݡ�
% 2.����������ֳ��Ⱥ� �����ʡ�
% 3.��������� (SNR) ���������
% 4.���ɽ�֯��������
% 5.Turbo �������ݡ�
% 6.Ӧ�� 16-QAM ���ƣ�����ƽ���źŹ��ʽ��й�һ����
% 7.�������ź�ͨ�� AWGN ͨ����
% 8.ʹ��LLR������������źţ��������أ�����ƽ���źŹ��ʽ��й�һ����
% 9.Turbo �������ݡ��������Խ������λӳ��˳���� Turbo ������Ԥ�ڵ�ӳ��˳���෴��
% 10.�������ͳ������
for k = 1:length(EbNo)
  errorStats = zeros(1,3); 
  noiseVar = 10^(-EbNo(k)/10)*(1/log2(M));
  chan.EbNo = EbNo(k);
  while errorStats(2) < 100 && errorStats(3) < 1e7
      data = randi([0 1],frmLen,1); % �����������������
      intrlvrInd = randperm(frmLen); % ��֯
      %����step������ִ�ж���
      encodedData = step(enc,data,intrlvrInd); % Turbo����
      modSignal = step(mod,encodedData); % ����
      receivedSignal = step(chan,modSignal); % AWGN�ŵ�
      demodSignal = step(demod,receivedSignal,noiseVar); % ���
      receivedBits = step(dec,-demodSignal,intrlvrInd); % ����
      errorStats = step(errorRate,data,receivedBits); % ͳ�Ʋ��
  end
  %����BER���ݲ����������ʶ���
  ber(k) = errorStats(1);
  reset(errorRate)
end
semilogy(EbNo,ber,'r-o')
grid %�л��ı��������ߵĿɼ���
xlabel('Eb/No (dB)')
ylabel('�������')
uncodedBER = berawgn(EbNo,'qam',M); % δ�����BER
hold on %����ԭͼ
semilogy(EbNo,uncodedBER,'b-')
legend('Turbo����','δ����')%ͼ��