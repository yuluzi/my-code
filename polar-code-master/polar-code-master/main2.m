%polar encoder and SC, BP , SCAN decoder for awgn channel

clear; tic;
global PCparams;


addpath('function');

N = 256;%码长
K = 8;%信息位长度
Rc = K/N;%码率
Rm = 1;%BPSK
ebn0 = 0:2:8;%;%一组信噪比0

design_snr_dB = 0;%使用BA构造方法构造极坐标代码的参数
sigma = 0.9;%使用GA构造方法构造极坐标代码的参数

min_frame_errors = 10;%30 最小错误帧数

max_frame_num = 100;%最大帧数
total_frame_num = 1e7;%任务总帧数

construction_method = input('choose polar code construction method (need run the file in constructed fold firstly to constructed the code): 0---BhattaBound  1---GA\n');
    % decoding_method = input('choose the decoding method: 1---SC  2---BP  3---SCAN  4---SCL 5---SSC\n');



ExecutorConditions = zeros(5,2);%用于存储各个执行体当前运行状态 0/1  0：不在运行 1：在运行
%第一列当前运行状态  其余列暂为用上

ExecutorNames = string(zeros(5));%存取每个执行体名称

ExecutorsBER = zeros(5,5);%存取各个执行体在不同的信噪比条件下执行的误码数
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码数

ExecutorsFER = zeros(5,5);%存取各个执行体在不同的信噪比条件下执行的误帧数
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧数

BER = zeros(5,5);%存取各个执行体在不同的信噪比条件下执行的误码率
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率

FER = zeros(5,5);%存取各个执行体在不同的信噪比条件下执行的误帧率
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧率

max_error_num = 10;%被裁决为错误执行体的最大次数


n=PCparams.n; 
disp('请输入各个译码器参数：');
bp_iter_num = input('\n input the iternum of the BP: at least 40 \n');
scan_iter_num = input('\n input the iternum of the SCAN: at least 1 \n');
scl_list_size = input('\n input the list size of the SCL: at least 1 \n');
crc_size = input('\n input the crc size of the SCL: at least 0 \n');
node_num = 2^(n+1)-1;
decoder_tree_initial = cell(1,node_num);
G_set = cell(1,n+1);
B_set = cell(1,n+1);


% 选取3个执行体构成执行体集合
executornum = 3; %input('请输入需要选取的执行体数量（请请输入最小为1，最大为5的正整数）：');
[executors,ExecutorConditions] = SwichExecutor(ExecutorConditions,executornum);   


for i = 1:length(executors)
    switch executors(i)
        case 1
            ExecutorNames(i) = ('1.SC decoder');
        case 2
            ExecutorNames(i) =('2.BP decoder');
        case 3
            ExecutorNames(i) =('3.SCAN decoder');
        case 4
            ExecutorNames(i) =('4.SCL decoder');
        case 5
            ExecutorNames(i) =('5.SSC decoder');
        otherwise
            break;
    end
     
end

%每个帧/块需要执行体集合中的每个执行体都执行一次
for ii = 1:length(executors)
    

    switch executors(ii)
        case 1
            fprintf('\n running polar_SC_decode now！');
            crc_size = 0;
            initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
        case 2
            fprintf('\n running polar_BP_decode now！');
            crc_size = 0;
            initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
        case 3
            fprintf('\n running polar_SCAN_decode now！');
            crc_size = 0;
            initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
        case 4
            fprintf('\nrunning polar_SCL_decode now！');
            initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
        case 5
            fprintf('\n running polar_SSC_decode now！');
            crc_size = 0;
            initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
            [decoder_tree_initial, G_set, B_set] = intial_tree_G( );
        otherwise
            fprintf('\n invalid input!!!');
            bp_iter_num = 60;
            scan_iter_num = 8;
            scl_list_size = 4;
            crc_size = 0;
    end


    %不同信噪比情况下分别执行
    for j = 1:length(ebn0)
            fprintf('\n Now running:%f  [%d of %d] \n\t Iteration-Counter: %53d',ebn0(j),j,length(ebn0),0);
            tt=tic();%用于启动一个计时器以测量代码的执行时间。返回值tt是一个计时器对象，用于停止计时器并计算经过的时间 
        %每次执行一个帧/块
        for l = 1:max_frame_num

            %运行 得到误码数和误帧数
            [ExecutorsBER(executors(ii),j),ExecutorsFER(executors(ii),j)] = RunABlock(N,K,executors(ii),ebn0(j),ExecutorsBER(executors(ii),j),ExecutorsFER(executors(ii),j),bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
            %计算误码率和误帧率


             fprintf(' \n异构执行体：%7s   在信噪比为：%7d  的环境下执行到%7d 个块的总误码数为：%7d ',ExecutorNames(ii),ebn0(j),l,ExecutorsBER(executors(ii),j));
             fprintf(' \n异构执行体：%7s   在信噪比为：%7d  的环境下执行到%7d 个块的总误帧数为：%7d ',ExecutorNames(ii),ebn0(j),l,ExecutorsFER(executors(ii),j));


            fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
            fprintf('\n');
        end
        
        %计算每个执行体在不同误码率下的误码率与误帧率
        BER(executors(ii),j) = ExecutorsBER(executors(ii),j)/(K*l);
        FER(executors(ii),j) = ExecutorsFER(executors(ii),j)/l;

    end


    % 画图

    subplot(2,1,1);%length(executors)个图片摆成一列，并放在做边一列
    xlabel('Eb/No (dB)')
    ylabel('polar decode BER');
    
    hold on; 
    switch executors(ii)
        case 1
            semilogy(ebn0,BER(executors(ii),1:5),'ro-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 2
            semilogy(ebn0,BER(executors(ii),1:5),'g*-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 3
            semilogy(ebn0,BER(executors(ii),1:5),'bd-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 4
            semilogy(ebn0,BER(executors(ii),1:5),'ys-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 5
            semilogy(ebn0,BER(executors(ii),1:5),'cp-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        otherwise
            fprintf('\n invalid input!!!');
    end
    hold on;

    subplot(2,1,2);%length(executors)个图片摆成两列，并放在右边一列
    %semilogy(ebn0,FER(executors(ii),1:5),'bs-','LineWidth',1.5,'MarkerSize',6);
    xlabel('Eb/No (dB)')
    ylabel('polar decode FER');
    
    hold on; 
    switch executors(ii)
        case 1
            semilogy(ebn0,FER(executors(ii),1:5),'ro-','LineWidth',1.5,'MarkerSize',6);
            hold on;  
        case 2
            semilogy(ebn0,FER(executors(ii),1:5),'g*-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 3
            semilogy(ebn0,FER(executors(ii),1:5),'bd-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 4
            semilogy(ebn0,FER(executors(ii),1:5),'ys-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        case 5
            semilogy(ebn0,FER(executors(ii),1:5),'cp-','LineWidth',1.5,'MarkerSize',6);
            hold on;
        otherwise
            fprintf('\n invalid input!!!');
    end
    hold on;
    subplot(2,1,1);
    switch executornum
        case 1
            legend(ExecutorNames(1));
            hold on;  
        case 2
            legend(ExecutorNames(1),ExecutorNames(2));
            hold on;
        case 3
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3));
            hold on;
        case 4
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3),ExecutorNames(4));
            hold on;
        case 5
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3),ExecutorNames(4),ExecutorNames(5));
            hold on;
    end
    grid on
    hold on;
    subplot(2,1,2);
    switch executornum
        case 1
            legend(ExecutorNames(1));
            hold on;  
        case 2
            legend(ExecutorNames(1),ExecutorNames(2));
            hold on;
        case 3
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3));
            hold on;
        case 4
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3),ExecutorNames(4));
            hold on;
        case 5
            legend(ExecutorNames(1),ExecutorNames(2),ExecutorNames(3),ExecutorNames(4),ExecutorNames(5));
            hold on;
    end
    grid on
    hold on;

end

toc

