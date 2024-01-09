% 执行体按照执行最大帧数出错概率为error_factor发生随机故障

clear; tic;

addpath('function');
global PCparams;
    
n=PCparams.n; 


% % 系统信噪比
% EbN0db = input(' 选择系统信噪比 Eb/N0（dB）,默认[0.5,1.0,1.5,2.0,2.5,3.0]: ');
% if isempty(EbN0db)
% EbN0db = [0.5,1.0,1.5,2.0,2.5,3.0]; % 默认
% end

N = 256;%码长
K = [64,128,192,224];%信息位长度

Rm = 1;%BPSK
ebn0 = 0:0.25:4;%;%一组信噪比0

%min_frame_error's = 10;%30 最小错误帧数
ExecutorNames = ["1.SC decoder","2.BP decoder","3.SCAN decoder","4.SCL decoder","5.SSC decoder"];%存取每个执行体名称
construction_method = input('请选择码构造方式: 0---BhattaBound  1---GA\n');
%decoding_method = input('choose the decoding method: 1---SC  2---BP  3---SCAN  4---SCL 5---SSC\n');
attack_method = input('请选择攻击方式: 0---无攻击  1---随机攻击  2---连续攻击 \n');
if attack_method == 2
    %选择一个被连续攻击的执行体
    rand_executor = randperm(4,1);
    fprintf(' \n异构执行体：%7s  将受到连续攻击\n ', ExecutorNames(rand_executor));
end

max_frame_num = input('请输入每个执行体集执行的最大帧数：\n');%最大帧数
total_frame_num = input('请输入本次任务共需执行的总帧数：\n');%任务总帧数
redundancy =  input('请输入DHR系统的冗余度（3 4 5）：\n');

SystemBERNum = 0; %系统误码总数，初始为0
SystemFERNum = 0; %系统误帧总数，初始为0

%初始化各个译码算法参数
disp('请输入各个译码器参数：');
bp_iter_num = input('\n input the iternum of the BP: at least 40 \n');
scan_iter_num = input('\n input the iternum of the SCAN: at least 1 \n');
scl_list_size = input('\n input the list size of the SCL: at least 1 \n');
crc_size = input('\n input the crc size of the SCL: at least 0 \n');
% node_num = 2^(n+1)-1;
% decoder_tree_initial = cell(1,node_num);
% G_set = cell(1,n+1);
% B_set = cell(1,n+1);
[decoder_tree_initial, G_set, B_set] = intial_tree_G( );

for jj = 1:length(K)
    Rc = K(jj)/N;%码率

    SystemBERNums = zeros(3,length(ebn0));%存取系统在不同的信噪比条件下执行的误码数
    % 例如：SystemBERNums(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误码数
    
    SystemFERNums = zeros(3,length(ebn0));%存取系统在不同的信噪比条件下执行的误帧数
    % 例如：SystemFERNums(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误帧数、
    
    SystemBERs = zeros(3,length(ebn0));%存取系统在不同的信噪比条件下执行的误码率
    % 例如：SystemBERs(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误码率
    
    SystemFERs = zeros(3,length(ebn0));%存取系统在不同的信噪比条件下执行的误帧率
    % 例如：SystemFERs(1,ebn0(1)) 代表系统在第一个信噪比条件下执行后得到的误帧率
     
    ComparativeExecutorsBERNum = zeros(5,length(ebn0)); %对比试验的各个执行体在不同的信噪比条件下执行的误码数量
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误码数
    
    ComparativeExecutorsFERNum = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误帧数量
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误帧数
    
    ComparativeExecutorsBER = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误码率
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误码率
    
    ComparativeExecutorsFER = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误帧率
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误帧率

    NoneAttackExecutorsBERNum = zeros(5,length(ebn0)); %对比试验的各个执行体在不同的信噪比条件下执行的误码数量
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误码数
    
    NoneAttackExecutorsFERNum = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误帧数量
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误帧数
    
    NoneAttackExecutorsBER = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误码率
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误码率
    
    NoneAttackExecutorsFER = zeros(5,length(ebn0));%对比试验的各个执行体在不同的信噪比条件下执行的误帧率
    % 例如：ComparativeExecutorsBERNum(executors(1),ebn0(1)) 代表对比试验的第一个执行体在第一个信噪比条件下执行后得到的误帧率
    
    SC_decoder_conditions = zeros(length(ebn0),total_frame_num);
    BP_decoder_conditions = zeros(length(ebn0),total_frame_num);
    SCAN_decoder_conditions = zeros(length(ebn0),total_frame_num);
    SCL_decoder_conditions = zeros(length(ebn0),total_frame_num);
    SSC_decoder_conditions = zeros(length(ebn0),total_frame_num);

    for j = 1:length(ebn0)
        switch attack_method
            case 0
                % DHR polar码仿真 
                 [ThreeSystemBERNum , ThreeSystemFERNum ,FourSystemBERNum , FourSystemFERNum,FiveSystemBERNum , FiveSystemFERNum, ...
                    ThreeSystemBER ,ThreeSystemFER, FourSystemBER ,FourSystemFER, FiveSystemBER,FiveSystemFER, TotalExecutorsBERNum , TotalExecutorsFERNum ,...
                    TotalExecutorsBER , TotalExecutorsFER , TotalNoneAttackExecutorsBERNum,TotalNoneAttackExecutorsFERNum,TotalNoneAttackExecutorsBER,...
                    TotalNoneAttackExecutorsFER,EachFrameExecutorsBERs ,EachFrameNoneAttackExecutorsBERs]   =     ...
                    NewSingleSNRSimulationRandomAttack(N,redundancy,K(jj),ebn0(j),total_frame_num,max_frame_num,construction_method,attack_method,   ...
                    bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set);
            case 1
                % DHR polar码仿真 
                [ThreeSystemBERNum , ThreeSystemFERNum ,FourSystemBERNum , FourSystemFERNum,FiveSystemBERNum , FiveSystemFERNum, ...
                    ThreeSystemBER ,ThreeSystemFER, FourSystemBER ,FourSystemFER, FiveSystemBER,FiveSystemFER, TotalExecutorsBERNum , TotalExecutorsFERNum ,...
                    TotalExecutorsBER , TotalExecutorsFER , TotalNoneAttackExecutorsBERNum,TotalNoneAttackExecutorsFERNum,TotalNoneAttackExecutorsBER,...
                    TotalNoneAttackExecutorsFER,EachFrameExecutorsBERs ,EachFrameNoneAttackExecutorsBERs]   =     ...
                    NewSingleSNRSimulationRandomAttack(N,redundancy,K(jj),ebn0(j),total_frame_num,max_frame_num,construction_method,attack_method,     ...
                    bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set);
            case 2
                % DHR polar码仿真 
                 [ThreeSystemBERNum , ThreeSystemFERNum ,FourSystemBERNum , FourSystemFERNum,FiveSystemBERNum , FiveSystemFERNum, ...
                    ThreeSystemBER ,ThreeSystemFER, FourSystemBER ,FourSystemFER, FiveSystemBER,FiveSystemFER, TotalExecutorsBERNum , TotalExecutorsFERNum ,...
                    TotalExecutorsBER , TotalExecutorsFER ,  EachFrameExecutorsBERs ,EachFrameNoneAttackExecutorsBERs]  =     ...
                    NewSingleSNRSimulationResistantAttack(N,redundancy,K(jj),ebn0(j),total_frame_num,max_frame_num,construction_method,     ...
                    bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set);
            otherwise
                fprintf('\n invalid input!!!');
        end
        
%          [SystemBERNum , SystemFERNum , SystemBER , SystemFER , TotalExecutorsBERNum , TotalExecutorsFERNum ,...
%             TotalExecutorsBER , TotalExecutorsFER ,  EachFrameExecutorsBERs]  =     ...
%             SingleSNRSimulationAttack2(N,K(jj),ebn0(j),total_frame_num,max_frame_num,construction_method,   ...
%             bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set);
%         switch redundancy
%             case 3
                SystemBERNums(1,j) = ThreeSystemBERNum;
                SystemFERNums(1,j) = ThreeSystemFERNum;
                SystemBERs(1,j) = ThreeSystemBER;
                SystemFERs(1,j) = ThreeSystemFER;
%             case 4
                SystemBERNums(2,j) = FourSystemBERNum;
                SystemFERNums(2,j) = FourSystemFERNum;
                SystemBERs(2,j) = FourSystemBER;
                SystemFERs(2,j) = FourSystemFER;
%             case 5
                SystemBERNums(3,j) = FiveSystemBERNum;
                SystemFERNums(3,j) = FiveSystemFERNum;
                SystemBERs(3,j) = FiveSystemBER;
                SystemFERs(3,j) = FiveSystemFER;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end

    
        fprintf(' \n 3余度的DHR系统本次运行共执行%7d 个块/帧，即%7d 个码数。在码率为：%7d ，信噪比为：%7d 的环境下 总误码数为：%7d ，总误帧数为：%7d ，系统误码率为：%7d ，系统误帧率为：%7d 。\n ',...
            total_frame_num , total_frame_num*K(jj) , Rc , ebn0(j) , SystemBERNums(1,j) , SystemFERNums(1,j) , SystemBERs(1,j) , SystemFERs(1,j));
    
        for i = 1:5
            ComparativeExecutorsBER(i,j) = TotalExecutorsBER(i,1);
            ComparativeExecutorsFER(i,j) = TotalExecutorsFER(i,1);
            ComparativeExecutorsBERNum(i,j) = TotalExecutorsBERNum(i,1);
            ComparativeExecutorsFERNum(i,j) = TotalExecutorsFERNum(i,1);
        end

        for i = 1:5
            NoneAttackExecutorsBER(i,j) = TotalNoneAttackExecutorsBER(i,1);
            NoneAttackExecutorsFER(i,j) = TotalNoneAttackExecutorsFER(i,1);
            NoneAttackExecutorsBERNum(i,j) = TotalNoneAttackExecutorsBERNum(i,1);
            NoneAttackExecutorsFERNum(i,j) = TotalNoneAttackExecutorsFERNum(i,1);
        end



        
        if K(jj) == 128 && ebn0(j) == 2
            figure(4);
            plot(1:1:total_frame_num,EachFrameExecutorsBERs(6,:),'k-','LineWidth',1.5,'MarkerSize',4)
            hold on;
            plot(1:1:total_frame_num,EachFrameExecutorsBERs(7,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
            hold on;
            plot(1:1:total_frame_num,EachFrameExecutorsBERs(8,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
            hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
            title('码率为1/2，信噪比为2时，各冗余度DHR系统的实时误码率');
            xlabel('时长（一帧运行的时间）')
            ylabel('误码率')
            legend('3-redundancy DHR system','4-redundancy DHR system','5-redundancy DHR system' )
        end





%         if K(jj) == 128 && ebn0(j) == 2
%             figure(4);
%             switch redundancy
%                 case 3
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(6,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 case 4
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(7,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 case 5
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(8,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 otherwise
%                     fprintf('\n redundancy input error!!!');
%             end
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameExecutorsBERs(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
%             title('码率为1/2，信噪比为2时，各编译码方法的实时误码率');
%             xlabel('时长（一帧运行的时间）')
%             ylabel('误码率')
%             legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%                 'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         end

%         if K(jj) == 128 && ebn0(j) == 2
%             figure(8);
%             switch redundancy
%                 case 3
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(6,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 case 4
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(7,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 case 5
%                     plot(1:1:total_frame_num,EachFrameExecutorsBERs(8,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                     hold on;
%                 otherwise
%                     fprintf('\n redundancy input error!!!');
%             end
%             plot(1:1:total_frame_num,EachFrameNoneAttackExecutorsBERs(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameNoneAttackExecutorsBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameNoneAttackExecutorsBERs(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameNoneAttackExecutorsBERs(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%             hold on;
%             plot(1:1:total_frame_num,EachFrameNoneAttackExecutorsBERs(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%             hold on;
%             title('码率为1/2，信噪比为2时，各编译码方法的实时误码率');
%             xlabel('时长（一帧运行的时间）')
%             ylabel('误码率')
%             legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%                 'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         end

   end

    
    if K(jj) == 64 
        figure(1);
        plot(ebn0,SystemBERs(1,:),'k-','LineWidth',1.5,'MarkerSize',4)
        hold on;
        plot(ebn0,SystemBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
        hold on;
        plot(ebn0,SystemBERs(3,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
        hold on;
        title('码率为1/4时，各编译码方法在不同信噪比环境下的误码率');
        xlabel('信噪比 (dB)')
        ylabel('误码率')
%         ylim([-0.01 1.01])
        legend('3-redundancy DHR system','4-redundancy DHR system','5-redundancy DHR system' )
        hold on;
        grid on;
    end

    if K(jj) == 128 
        figure(2);
        plot(ebn0,SystemBERs(1,:),'k-','LineWidth',1.5,'MarkerSize',4)
        hold on;
        plot(ebn0,SystemBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
        hold on;
        plot(ebn0,SystemBERs(3,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
        hold on;
        title('码率为1/2时，各编译码方法在不同信噪比环境下的误码率');
        xlabel('信噪比 (dB)')
        ylabel('误码率')
%         ylim([-0.01 1.01])
        legend('3-redundancy DHR system','4-redundancy DHR system','5-redundancy DHR system')
        hold on;
        grid on;
    end

    if K(jj) == 192
        figure(3);
        plot(ebn0,SystemBERs(1,:),'k-','LineWidth',1.5,'MarkerSize',4)
        hold on;
        plot(ebn0,SystemBERs(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
        hold on;
        plot(ebn0,SystemBERs(3,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
        hold on;
        title('码率为3/4时，各编译码方法在不同信噪比环境下的误码率');
        xlabel('信噪比 (dB)')
        ylabel('误码率')
%         ylim([-0.01 1.01])
        legend('3-redundancy DHR system','4-redundancy DHR system','5-redundancy DHR system')
        hold on;
        grid on;
    end



















%     if K(jj) == 64 
%         figure(1);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,Syste-mBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,ComparativeExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为1/4时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end
% 
%     if K(jj) == 128 
%         figure(2);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,SystemBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,ComparativeExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为1/2时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end
% 
%     if K(jj) == 192
%         figure(3);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,SystemBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,ComparativeExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,ComparativeExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为3/4时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end
% 
% 
%     if K(jj) == 64 
%         figure(5);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,SystemBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,NoneAttackExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为1/4时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end
% 
%     if K(jj) == 128 
%         figure(6);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,SystemBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,NoneAttackExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为1/2时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end
% 
%     if K(jj) == 192
%         figure(7);
%         switch redundancy
%             case 3
%                 plot(ebn0,SystemBERs(1,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 4
%                 plot(ebn0,SystemBERs(2,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             case 5
%                 plot(ebn0,SystemBERs(3,:),'k-.','LineWidth',1.5,'MarkerSize',4)
%                 hold on;
%             otherwise
%                 fprintf('\n redundancy input error!!!');
%         end
%         plot(ebn0,NoneAttackExecutorsBER(1,:),'-','LineWidth',1.5,'Color',[0.25,0.25,0.25])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(2,:),'-.','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(3,:),'-','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(4,:),':','LineWidth',1.5,'Color',[0.5,0.5,0.5])
%         hold on;
%         plot(ebn0,NoneAttackExecutorsBER(5,:),'--','LineWidth',1.5,'Color',[0.75,0.75,0.75])
%         hold on;
%         title('码率为3/4时，各编译码方法在不同信噪比环境下的误码率');
%         xlabel('信噪比 (dB)')
%         ylabel('误码率')
% %         ylim([-0.01 1.01])
%         legend('BER of DHR polar decoder','BER of SC polar decoder','BER of BP polar decoder' ,...
%             'BER of SCAN polar decoder','BER of SCL polar decoder','BER of SSC polar decoder')
%         hold on;
%         grid on;
%     end

end

    %画图 

    



toc

