%polar encoder and SC, BP , SCAN decoder for awgn channel

clear; tic;
global PCparams;
addpath('function');


N = 256;%码长
K = 8;%信息位长度
Rc = K/N;%码率
Rm = 1;%BPSK
ebn0 = 0:2:8;%一组信噪比
ebn0_num = 10.^(ebn0/10);%与信噪比对应的信噪比值
SNR = ebn0 + 10*log10(Rc*Rm)+10*log10(2);%与信噪比对应的信道信噪比
noise_sigma = 1./(10.^(SNR/10));%信噪比对应的噪声方差 

design_snr_dB = 0;%使用BA构造方法构造极坐标代码的参数
sigma = 0.9;%使用GA构造方法构造极坐标代码的参数

min_frame_errors = 10;%30 最小错误帧数
min_frame_num = 1000;%最小帧数
max_frame_num = 1e7;%最大帧数

%datas = zeros(5,3);%用于存储各个异构体的BER、FER以及
%第一列是BER 第二列是FER 第三列是记录每个执行体每次执行的块错误率 第四列是

construction_method = input('choose polar code construction method (need run the file in constructed fold firstly to constructed the code): 0---BhattaBound  1---GA\n');
    % decoding_method = input('choose the decoding method: 1---SC  2---BP  3---SCAN  4---SCL 5---SSC\n');


executor =  randperm(5,1)-1;%产生2个1-4之间的随机数
disp('选取的异构执行体编号如下：');
disp(executor);%输出生成的随机数组

%循环执行各个执行体
for i = 1:length(executor)
    switch executor(i)
    case 1
        disp('running polar_SC_decode now！');
        crc_size = 0;
        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
    case 2
        disp('running polar_BP_decode now！');
        bp_iter_num = input('input the iternum of the BP: at least 40 \n');
        crc_size = 0;
        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
    case 3
        disp('running polar_SCAN_decode now！');
        scan_iter_num = input('input the iternum of the SCAN: at least 1 \n');
        crc_size = 0;
        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
    case 4
        disp('running polar_SCL_decode now！');
        scl_list_size = input('input the list size of the SCL: at least 1 \n');
        crc_size = input('input the crc size of the SCL: at least 0 \n');
        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
    case 5
        disp('running polar_SSC_decode now！');
        crc_size = 0;
        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
        [decoder_tree_initial, G_set, B_set] = intial_tree_G( );
    otherwise
        disp('invalid input!!!');
        bp_iter_num = 60;
        scan_iter_num = 8;
        scl_list_size = 4;
        crc_size = 0;
    end



    n = PCparams.n;
    F = [1 0;1 1];
    B=1;
    for ii=1:n
        B = kron(B,F);
    end
    F_kron_n = B;%克罗内克矩阵
    
    FER = zeros(1,length(ebn0));%初始化误码率
    BER = zeros(1,length(ebn0));%初始化比特错误率
    bpsk_FER=zeros(1,length(ebn0));%初始化bpsk误码率       
    bpsk_BER=zeros(1,length(ebn0));%初始化bpsk比特错误率   
    
    for j = 1:length(ebn0)
	    fprintf('\n Now running:%f  [%d of %d] \n\t Iteration-Counter: %53d',ebn0(j),j,length(ebn0),0);
	    tt=tic();%用于启动一个计时器以测量代码的执行时间。返回值tt是一个计时器对象，用于停止计时器并计算经过的时间
	    for l = 1:max_frame_num
		    
            de_bpsk = zeros(1,N);
            u=randi(2,1,K)-1; %Bernoulli(0.5);  伯努利分布 生成一个长度为K的随机二进制序列（0或1），每个元素的取值都是等概率的
		    x=pencode(u,F_kron_n);  %编码程序
            tx_waveform=bpsk(x); % bpsk调制
            noise=sqrt(noise_sigma(j))*randn(1,N); % awgn噪声
            rx_waveform = tx_waveform+noise;
            de_bpsk(rx_waveform>0)=1; % bpsk解调
                
            % nfails 表示传输过程中发生错误的比特数，即解调后得到的比特与原始发送的比特不相等的数量。
            % de_bpsk 表示解调后的比特序列，x 表示发送的比特序列。
            % 因此，sum(de_bpsk ~= x) 表示解调后比特序列中与原始发送的比特序列不相等的比特数之和，也就是错误的比特数。
            nfails = sum(de_bpsk ~= x);
            bpsk_FER(j) = bpsk_FER(j) + (nfails>0);
            bpsk_BER(j) = bpsk_BER(j) + nfails;
            
            initia_llr = -2*rx_waveform/noise_sigma(j);
	    
            switch executor(i)             
                case 1
                    [u_llr] = polar_SC_decode(initia_llr);
                case 2
                    [u_llr,~] = polar_BP_decode(initia_llr,bp_iter_num);
                case 3
                    [u_llr,~] = polar_SCAN_decode(initia_llr,scan_iter_num);
                case 4
                    [u_llr] = polar_SCL_decode(initia_llr,scl_list_size);      
                case 5
                    u_hard_decision = polar_SSC_decode(decoder_tree_initial, G_set, B_set, initia_llr);
                    u_llr = 1-2*u_hard_decision;
            end
            if PCparams.crc_size
                uhat_crc_llr = u_llr(PCparams.FZlookup == -1)';
                uhat_llr = uhat_crc_llr (1:PCparams.K);
            else
                uhat_llr = u_llr(PCparams.FZlookup == -1)';
            end
		    uhat = zeros(1,K);
            uhat(uhat_llr<0) =1;
    
		    nfails = sum(uhat ~= u);
            FER(j) = FER(j) + (nfails>0);
            BER(j) = BER(j) + nfails;
		    if mod(l,20)==0
                for iiiii=1:53
                    fprintf('\b');
                end
                fprintf(' %7d   ---- %7d FEs, %7d BEs found so far',l,FER(j),BER(j));
            end
		    if l>=min_frame_num && FER(j)>=min_frame_errors  %frame errors, sufficient to stop
                break;
            end
	    end
            
        BER(j) = BER(j)/(K*l);
        FER(j) = FER(j)/l;

        %将各个执行的BER和FER存储到datas数组中
        %datas(executor(i),1) = BER(j);
        %datas(executor(i),1) = FER(j);
        
        bpsk_BER(j) = bpsk_BER(j)/(N*l);
        bpsk_FER(j) = bpsk_FER(j)/l;
    
        fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
        fprintf('\n');
    end

    subplot(length(executor),1,i);%length(executor)个图片摆成一列
    semilogy(ebn0,BER,'bs-','LineWidth',1.5,'MarkerSize',6)
    xlabel('Eb/No (dB)')
    ylabel('BER/FER')
    hold on;
    
    semilogy(ebn0,FER,'gs:','LineWidth',1.5,'MarkerSize',6)
    hold on;
    
    semilogy(ebn0,bpsk_BER,'rv-','LineWidth',1.5,'MarkerSize',6)
    hold on;
    
    semilogy(ebn0,bpsk_FER,'rv:','LineWidth',1.5,'MarkerSize',6)
    hold on;
    
    grid on;
    %legend('PC SCL-8-16 BER','PC SCL-8-16 FER','BPSK BER','BPSK FER');


    switch executor(i)
        case 1
            %legend('polar SC decode BER','BPSK BER');
            legend('polar SC decode BER','polar SC decode FER','BPSK BER','BPSK FER');
        case 2
            %legend('polar BP decode BER','BPSK BER');
            legend('polar BP decode BER','polar BP decode FER','BPSK BER','BPSK FER');
        case 3
             %legend('polar SCAN decode BER','BPSK BER');
            legend('polar SCAN decode BER','polar SCAN decode FER','BPSK BER','BPSK FER');
        case 4
            %legend('polar SCL decode BER','BPSK BER');
            legend('polar SCL decode BER','polar SCL decode FER','BPSK BER','BPSK FER');
        case 5
            %legend('polar SSC decode BER','BPSK BER');
            legend('polar SSC decode BER','polar SSC decode FER','BPSK BER','BPSK FER');
    end
    

end



toc

