%polar encoder and SC, BP , SCAN decoder for awgn channel

clear; tic;
global PCparams;
addpath('function');


N = 256;%码长
K = 8;%信息位长度
Rc = K/N;%码率
Rm = 1;%BPSK


design_snr_dB = 0;%使用BA构造方法构造极坐标代码的参数
sigma = 0.9;%使用GA构造方法构造极坐标代码的参数

min_frame_errors = 10;%30 最小错误帧数
min_frame_num = 10000;%最小帧数
max_frame_num = 1e7;%最大帧数

ExecutorConditions = zeros(5,5);%用于存储各个执行体被判决成错误执行体的次数 以及 当前运行状态 0/1  0：不在运行 1：在运行
%第一列是记录每个执行体每次执行的块错误率 第二列是当前运行状态 

ExecutorsBER = zeros(5,3);%存取各个执行体在不同的信噪比条件下执行的误码率 
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率

ExecutorsFER = zeros(5,5);%存取各个执行体在不同的信噪比条件下执行的误帧率 
% 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率

max_error_num = 10;%被裁决为错误执行体的最大次数
construction_method = input('choose polar code construction method (need run the file in constructed fold firstly to constructed the code): 0---BhattaBound  1---GA\n');
    % decoding_method = input('choose the decoding method: 1---SC  2---BP  3---SCAN  4---SCL 5---SSC\n');

 
% n = PCparams.n;
% F = [1 0;1 1];
% B=1;
% for ii=1:n
%     B = kron(B,F);
% end
% F_kron_n = B;%克罗内克矩阵

% FER = zeros(1,length(ebn0));%初始化误码率
% BER = zeros(1,length(ebn0));%初始化比特错误率
% bpsk_FER=zeros(1,length(ebn0));%初始化bpsk误码率
% bpsk_BER=zeros(1,length(ebn0));%初始化bpsk比特错误率

% 选取三个执行体构成执行体集合
executors = [SwichExecutor(ExecutorConditions),SwichExecutor(ExecutorConditions),SwichExecutor(ExecutorConditions)];

%不同信噪比情况下分别执行
% for j = 1:length(ebn0)
%     fprintf('\n Now running:%f  [%d of %d] \n\t Iteration-Counter: %53d',ebn0(j),j,length(ebn0),0);
%     tt=tic();%用于启动一个计时器以测量代码的执行时间。返回值tt是一个计时器对象，用于停止计时器并计算经过的时间

    %每次执行一个帧/块
    for l = 1:max_frame_num

        %每个帧/块需要执行体集合中的每个执行体都执行一次
        for ii = 1:length(executors)
            
            RunABlock(N,K,construction_method,design_snr_dB,sigma,executors(ii),ExecutorsBER,ExecutorsFER);
            
            
%             %将各个执行的BER和FER存储到datas数组中
%             ExecutorConditions(executor(ii),1) = BER(j);
%             ExecutorConditions(executor(ii),1) = FER(j);
            
            fprintf(' \n异构执行体：%7d   在信噪比为：%7d  的环境下执行后的误码率为：%7d ',executors(ii),ebn0(j),ExecutorsBER(executors(ii),j));
            fprintf(' \n异构执行体：%7d   在信噪比为：%7d  的环境下执行后的误帧率为：%7d ',executors(ii),ebn0(j),ExecutorsFER(executors(ii),j));
            %disp('异构执行体：',executors(ii),'在信噪比为：',ebn0(j),'的环境下执行后的误码率为：',ExecutorsBER(executors(ii),j));
            %disp('异构执行体：',executors(ii),'在信噪比为：',ebn0(j),'的环境下执行后的误帧率为：',ExecutorsBER(executors(ii),j));

        
        end




%         if mod(l,20)==0
%             for iiiii=1:53
%                 fprintf('\b');
%             end
%             fprintf(' %7d   ---- %7d FEs, %7d BEs found so far',l,FER(j),BER(j));
%         end
%         if l>=min_frame_num && FER(j)>=min_frame_errors  %frame errors, sufficient to stop
%             break;
%         end


%         BER = BER(j)/(K*l);
%         FER = FER(j)/l;

    end

    for i = 1:5
        for j = 1:2
            fprintf(' \n异构执行体：%7d   在信噪比为：%7d  的环境下执行后的误码率为：%7d ',i,ebn0(j),ExecutorsBER(i,j));
        end
    end

%

%     fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
%     fprintf('\n');
% end


%画图
% subplot(length(executor),1,i);%length(executor)个图片摆成一列
% semilogy(ebn0,BER,'bs-','LineWidth',1.5,'MarkerSize',6)
% xlabel('Eb/No (dB)')
% ylabel('BER/FER')
% hold on;
% 
% semilogy(ebn0,FER,'gs:','LineWidth',1.5,'MarkerSize',6)
% hold on;
% 
% semilogy(ebn0,bpsk_BER,'rv-','LineWidth',1.5,'MarkerSize',6)
% hold on;
% 
% semilogy(ebn0,bpsk_FER,'rv:','LineWidth',1.5,'MarkerSize',6)
% hold on;
% 
% grid on;
% %legend('PC SCL-8-16 BER','PC SCL-8-16 FER','BPSK BER','BPSK FER');
% 
% 
% switch executor(i)
%     case 1
%         %legend('polar SC decode BER','BPSK BER');
%         legend('polar SC decode BER','polar SC decode FER','BPSK BER','BPSK FER');
%     case 2
%         %legend('polar BP decode BER','BPSK BER');
%         legend('polar BP decode BER','polar BP decode FER','BPSK BER','BPSK FER');
%     case 3
%         %legend('polar SCAN decode BER','BPSK BER');
%         legend('polar SCAN decode BER','polar SCAN decode FER','BPSK BER','BPSK FER');
%     case 4
%         %legend('polar SCL decode BER','BPSK BER');
%         legend('polar SCL decode BER','polar SCL decode FER','BPSK BER','BPSK FER');
%     case 5
%         %legend('polar SSC decode BER','BPSK BER');
%         legend('polar SSC decode BER','polar SSC decode FER','BPSK BER','BPSK FER');
% end



toc

