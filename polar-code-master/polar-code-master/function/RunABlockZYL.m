function  [BER,FER] = RunABlockZYL(N,K,u,executor,ebn0,ThisExecutorBER,ThisExecutorFER,...
    bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set,attack_t,ran_t,rand_executoers)
% 本函数用于模拟一个执行体在某一个特定误码率环境下执行一个帧的过程，并更新本次执行后的误码数与误帧数累加和
%
% 输入：
% N:码长
% K:信息位长度
% executor:当前执行的异构执行体编号 
% ebn0：当前执行时的信噪比 
% ThisExecutorBER：当前执行体在该信噪比环境下执行得到的误码数累加和 
% ThisExecutorFER：当前执行体在该信噪比环境下执行得到的误帧数累加和 
% bp_iter_num：BP译码参数 
% scan_iter_num：SCAN译码参数
% scl_list_size：SCL译码参数
% decoder_tree_initial：SSC译码参数，用于生成译码树结构
% G_set：SSC译码参数，用于初始化SSC译码时将用到的G矩阵
% B_set：SSC译码参数，用于初始化SSC译码时将用到的B_N矩阵
%
% 输出：
% BER：当前执行体在该信噪比环境下执行得到的误码数累加和，即ThisExecutorBER中的值加上这一个帧的误码数
% FER：当前执行体在该信噪比环境下执行得到的误帧数累加和，即ThisExecutorFER中的值加上这一个帧的误帧数

    
    %定义BER和FER的初始值,分别代表这一个执行体在某一个信噪比条件下运行到第 l 个帧时的误码数和误帧数
    BER = ThisExecutorBER;
    FER = ThisExecutorFER;


    Rc = K/N;%码率
    Rm = 1;%BPSK
    %ebn0 = 0:2:8;%一组信噪比
    %ebn0_num = 10.^(ebn0/10);%与信噪比对应的信噪比值
    SNR = ebn0 + 10*log10(Rc*Rm)+10*log10(2);%与信噪比对应的信道信噪比
    noise_sigma = 1./(10.^(SNR/10));%信噪比对应的噪声方差 

    global PCparams;
    n = PCparams.n;
    F = [1 0;1 1];
    B=1;
    for ii=1:n
        B = kron(B,F);
    end
    F_kron_n = B;%克罗内克矩阵

        %de_bpsk = zeros(1,N);
        x=pencode(u,F_kron_n);  %编码程序

        %加人为干扰（随机选个时间段挨打，找两个执行体）
        if executor==rand_executoers && ran_t >= attack_t 
           num_errors = ceil(N*0.1);%10%的误码率
           error_indices = randperm(N,num_errors);
           for kk = 1:num_errors
               x(error_indices(kk)) = mod(x(error_indices(kk))+1,2);
           end
        end
       
        tx_waveform=bpsk(x); % bpsk调制
        noise=sqrt(noise_sigma)*randn(1,N); % awgn噪声
        rx_waveform = tx_waveform+noise;
        %de_bpsk(rx_waveform>0)=1; % bpsk解调
        
        % nfails 表示传输过程中发生错误的比特数，即解调后得到的比特与原始发送的比特不相等的数量。
        % de_bpsk 表示解调后的比特序列，x 表示发送的比特序列。
        % 因此，sum(de_bpsk ~= x) 表示解调后比特序列中与原始发送的比特序列不相等的比特数之和，也就是错误的比特数。
        %nfails = sum(de_bpsk ~= x); %bpsk的nfails
        
        % bpsk_FER(j) = bpsk_FER(j) + (nfails>0);
        % bpsk_BER(j) = bpsk_BER(j) + nfails;
        
        
        initia_llr = -2*rx_waveform/noise_sigma;
        
        switch executor
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

        FER = FER + (nfails>0);
        BER = BER + nfails;


end

