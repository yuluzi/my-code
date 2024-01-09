function [SystemBERs,SystemFERs] = mainZYL()
%polar encoder and SC, BP , SCAN decoder for awgn channel
% 每个执行体集合在一个信噪比环境下执行一个max_frame_num就重新选取执行体
clear; tic;

global PCparams;
n = PCparams.n;
addpath('function');

N = 256;%码长
K =128;%信息位长度
Rc = K/N;%码率
Rm = 1;%BPSK
ebn0 = 0:0.5:4.5;%;%一组信噪比0




%min_frame_errors = 10;%30 最小错误帧数

construction_method = input('choose polar code construction method (need run the file in constructed fold firstly to constructed the code): 0---BhattaBound  1---GA\n');
%decoding_method = input('choose the decoding method: 1---SC  2---BP  3---SCAN  4---SCL 5---SSC\n');

max_frame_num = input('请输入每个执行体集执行的最大帧数：');%最大帧数
total_frame_num = input('请输入本次任务共需执行的总帧数：');%任务总帧数

SystemBERNum = 0; %系统误码总数，初始为0
SystemFERNum = 0;   %系统误帧总数，初始为0

%初始化各个译码算法参数
disp('请输入各个译码器参数：');
bp_iter_num = input('\n input the iternum of the BP: at least 40 \n');
scan_iter_num = input('\n input the iternum of the SCAN: at least 1 \n');
scl_list_size = input('\n input the list size of the SCL: at least 1 \n');
crc_size = input('\n input the crc size of the SCL: at least 0 \n');
node_num = 2^(n+1)-1;
decoder_tree_initial = cell(1,node_num);
G_set = cell(1,n+1);
B_set = cell(1,n+1);

SystemBERNums = zeros(1,length(ebn0));%存取系统在不同的信噪比条件下执行的误码数
% 例如：SystemBERNums(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误码数

SystemFERNums = zeros(1,length(ebn0));%存取系统在不同的信噪比条件下执行的误帧数
% 例如：SystemFERNums(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误帧数、

SystemBERs = zeros(1,length(ebn0));%存取系统在不同的信噪比条件下执行的误码率
% 例如：SystemBERs(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误码率

SystemFERs = zeros(1,length(ebn0));%存取系统在不同的信噪比条件下执行的误帧率
% 例如：SystemFERs(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误帧率


for j = 1:length(ebn0)
    
    [SystemBERNum , SystemFERNum , SystemBER , SystemFER] =     ...
        SingleSNRSimulation(N,K,ebn0(j),total_frame_num,max_frame_num,construction_method,   ...
        bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set);

    SystemBERNums(1,j) = SystemBERNum;
    SystemFERNums(1,j) = SystemFERNum;
    SystemBERs(1,j) = SystemBER;
    SystemFERs(1,j) = SystemFER;

    fprintf(' \n 本次运行共执行%7d 个块/帧，即%7d 个码数。在信噪比为：%7d 的环境下 总误码数为：%7d ，总误帧数为：%7d ，系统误码率为：%7d ，系统误帧率为：%7d 。\n ',...
        total_frame_num , total_frame_num*K , ebn0(j) , SystemBERNums(1,j) , SystemFERNums(1,j) , SystemBERs(1,j) , SystemFERs(1,j));
end


% 画图
%subplot(2,1,1);     % 2个图片摆成1列，第三个参数是指第几个小图
%semilogy(ebn0,SystemBERs,'bs-','LineWidth',1.5,'MarkerSize',6)
%xlabel('Eb/No (dB)')
%ylabel('SystemBER')
%legend('BER of polar decode in different ebn0s')
%hold on;
%grid on;

%subplot(2,1,2);
%semilogy(ebn0,SystemFERs,'gs:','LineWidth',1.5,'MarkerSize',6)
%xlabel('Eb/No (dB)')
%ylabel('SystemFER')
%legend('FER of polar decode in different ebn0s')
%hold on;
%grid on;


toc

