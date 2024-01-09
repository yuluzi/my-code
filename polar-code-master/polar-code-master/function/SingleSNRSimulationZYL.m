function [SystemBERNum , SystemFERNum , SystemBER , SystemFER] =  ...
SingleSNRSimulation(N,K,ebn0,total_frame_num,max_frame_num,construction_method,   ...
bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set)
%SINGLESNRSIMULATION 此处显示有关此函数的摘要
%   此处显示详细说明



    
    SystemBERNum = 0;
    SystemFERNum = 0;
    design_snr_dB = 0;%使用BA构造方法构造极坐标代码的参数
    sigma = 0.9;%使用GA构造方法构造极坐标代码的参数


    
    IsContinue = true;
    current_total_frame_num = 0;
    
    ExecutorConditions = zeros(5,2);%用于存储各个执行体当前运行状态 0/1  0：不在运行 1：在运行
    %第一列当前运行状态  其余列暂未用上
    
    ExecutorNames = ["1.SC decoder","2.BP decoder","3.SCAN decoder","4.SCL decoder","5.SSC decoder"];%存取每个执行体名称
    %选择一个被攻击执行体
    rand_executor = randperm(5,1);
  
    while IsContinue
        if current_total_frame_num>=total_frame_num
            break;
        else
            ExecutorsBERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误码数
            % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码数
            
            ExecutorsFERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误帧数
            % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧数
            
            ExecutorsBER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误码率
            % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率
            
            ExecutorsFER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误帧率
            % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧率
            
            
            % 选取3个执行体构成执行体集合
            executornum = 3; %input('请输入需要选取的执行体数量（请请输入最小为1，最大为5的正整数）：');
            [executors,ExecutorConditions] = SwichExecutor(ExecutorConditions,executornum);   
            
            
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
                %生成随机攻击时间数组
                attack_ts = rand(1,10);
                %不同信噪比情况下分别执行
                for j = 1:length(ebn0)
                        fprintf('\n Now running:%f  [%d of %d] \n\t Iteration-Counter: %53d',ebn0(j),j,length(ebn0),0);
                        tt=tic();%用于启动一个计时器以测量代码的执行时间。返回值tt是一个计时器对象，用于停止计时器并计算经过的时间 
                    %每次执行一个指定帧/块数
                    for l = current_total_frame_num + 1 : current_total_frame_num + max_frame_num
                        %每次先生成信息序列，然后执行，使每个执行体执行相同的信息序列
                        u = randi(2,1,K)-1; %Bernoulli(0.5);  伯努利分布 生成一个长度为K的随机二进制序列（0或1），每个元素的取值都是等概率的
                        %运行 得到误码数和误帧数
                        [ExecutorsBERNum(executors(ii),j),ExecutorsFERNum(executors(ii),j)] = ...
                            RunABlock(N,K,u,executors(ii),ebn0(j),ExecutorsBERNum(executors(ii),j),ExecutorsFERNum(executors(ii),j),...
                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set,attack_ts(j),toc(tt),rand_executor);
                        %计算误码数和误帧数
                        if l ==current_total_frame_num + max_frame_num
                            fprintf(' \n异构执行体：%7s   在信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误码数为：%7d ',ExecutorNames(executors(ii)),ebn0(j),l-max_frame_num+1,l,ExecutorsBERNum(executors(ii),j));
                            fprintf(' \n异构执行体：%7s   在信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误帧数为：%7d ',ExecutorNames(executors(ii)),ebn0(j),l-max_frame_num+1,l,ExecutorsFERNum(executors(ii),j));
                            fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
                            fprintf('\n');
                        end
                        
                    end
                    
                    %计算每个执行体在不同信噪比下的误码率与误帧率
                    ExecutorsBER(executors(ii),j) = ExecutorsBER(executors(ii),j)/(K*l);
                    ExecutorsFER(executors(ii),j) = ExecutorsFER(executors(ii),j)/l;
    
                end
            
            end
    
            %所有执行体执行完一轮（max_frame_num）后，比较被选出的执行体的BER和FER，并选最好的作为系统的结果,最差的下次不被选用
            [MINBERs,MINBERExecutors] = GetMinBER(length(executors),ExecutorsBERNum,ebn0,executors,ExecutorNames);
            [MINFERs,MINFERExecutors] = GetMinFER(length(executors),ExecutorsFERNum,ebn0,executors,ExecutorNames);
            [MAXBERs,MAXBERExecutors] = GetMaxBER(length(executors),ExecutorsBERNum,ebn0,executors,ExecutorNames);
            [MAXFERs,MAXFERExecutors] = GetMaxFER(length(executors),ExecutorsFERNum,ebn0,executors,ExecutorNames);
            fprintf(' \n异构执行体：%7s   不会在下一轮执行体选取中被选中 ',ExecutorNames(MAXBERExecutors));
            
            %ExecutorConditions全部清零，以便下次选取
            for i = 1:5
                ExecutorConditions(i,1) = 0;
            end

            if MINBERExecutors ~= MAXBERExecutors %排除误码数为0时剔除一个异构体
                ExecutorConditions(MAXBERExecutors,1) = 1;
            end

            %计算系统的误码总数和误帧总数
            SystemBERNum = SystemBERNum + MINBERs;
            SystemFERNum = SystemFERNum + MINFERs;
            SystemBER = SystemBERNum/(total_frame_num*K);
            SystemFER = SystemFERNum/total_frame_num;
    
    
            current_total_frame_num = current_total_frame_num + max_frame_num;
    
        end
       
    end





end

