function [ThreeSystemBERNum , ThreeSystemFERNum ,FourSystemBERNum , FourSystemFERNum,FiveSystemBERNum , FiveSystemFERNum, ...
    ThreeSystemBER ,ThreeSystemFER, FourSystemBER ,FourSystemFER, FiveSystemBER,FiveSystemFER, TotalExecutorsBERNum , TotalExecutorsFERNum ,...
    TotalExecutorsBER , TotalExecutorsFER , TotalNoneAttackExecutorsBERNum,TotalNoneAttackExecutorsFERNum,TotalNoneAttackExecutorsBER,...
    TotalNoneAttackExecutorsFER,EachFrameExecutorsBERs ,EachFrameNoneAttackExecutorsBERs] =  ...
    NewSingleSNRSimulationRandomAttack(N,redundancy,K,ebn0,total_frame_num,max_frame_num,construction_method,attack_method,     ...
    bp_iter_num,scan_iter_num,scl_list_size,crc_size,decoder_tree_initial, G_set, B_set)
% SingleSNRSimulationDHR 此处显示有关此函数的摘要
% 思路：本函数每个执行体跑两次，一次被攻击，一次不被攻击，每次系统输出分别去相应的执行体中拉数据裁决
%
%   执行体按照执行最大帧数出错概率为error_factor发生随随机故障

% r = a + (b-a)*rand(m,n) = 0.3 + （1 - 0.3）* rand(1)

    ThreeSystemBERNum = 0;
    ThreeSystemFERNum = 0;

    FourSystemBERNum = 0;
    FourSystemFERNum = 0;

    FiveSystemBERNum = 0;
    FiveSystemFERNum = 0;
    
    ThreeSystemBER = 0;
    ThreeSystemFER = 0;

    FourSystemBER = 0;
    FourSystemFER = 0;

    FiveSystemBER = 0;
    FiveSystemFER = 0;

%     SystemBERNum = 0;
%     SystemFERNum = 0;
    design_snr_dB = 0;%使用BA构造方法构造极坐标代码的参数
    sigma = 0.9;%使用GA构造方法构造极坐标代码的参数
    
    switch attack_method
        case 0
            % 生成执行体攻击频率 保留小数点后一位
            SC_attack_rate = 0;%roundn(0.1 + (0.22 - 0.1)* rand(1) , -1);
            BP_attack_rate = 0;%roundn(0.1 + (0.35 - 0.1)* rand(1) , -1);
            SCAN_attack_rate = 0;%roundn(0.1 + (0.28 - 0.1)* rand(1) , -1);
            SCL_attack_rate = 0;%roundn(0.1 + (0.33 - 0.1)* rand(1) , -1);
            SSC_attack_rate = 0;%roundn(0.1 + (0.18 - 0.1)* rand(1) , -1);
        case 1
            % 生成执行体攻击频率 保留小数点后一位
            SC_attack_rate = roundn(0.1 + (0.22 - 0.1)* rand(1) , -1);
            BP_attack_rate = roundn(0.1 + (0.35 - 0.1)* rand(1) , -1);
            SCAN_attack_rate = roundn(0.1 + (0.28 - 0.1)* rand(1) , -1);
            SCL_attack_rate = roundn(0.1 + (0.33 - 0.1)* rand(1) , -1);
            SSC_attack_rate = roundn(0.1 + (0.18 - 0.1)* rand(1) , -1);
        otherwise
            fprintf('\n invalid input!!!');
                
    end
    
    IsContinue = true;
    current_total_frame_num = 0;
    
    ExecutorConditions = zeros(5,2);%用于存储各个执行体当前运行状态 -1/0/1  -1：故障   0：不在运行   1：在运行
    %第一列当前运行状态  其余列暂未用上
    SC_decoder_conditions = zeros(1,total_frame_num);
    BP_decoder_conditions = zeros(1,total_frame_num);
    SCAN_decoder_conditions = zeros(1,total_frame_num);
    SCL_decoder_conditions = zeros(1,total_frame_num);
    SSC_decoder_conditions = zeros(1,total_frame_num);


    ExecutorNames = ["1.SC decoder","2.BP decoder","3.SCAN decoder","4.SCL decoder","5.SSC decoder"];%存取每个执行体名称
    
    TotalExecutorsBERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误码数
    % 例如：TotalExecutorsBERNum(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码数
    
    TotalExecutorsFERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误帧数
    % 例如：TotalExecutorsFERNum(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧数
    
    TotalExecutorsBER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误码率
    % 例如：TotalExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率
    
    TotalExecutorsFER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误帧率
    % 例如：TotalExecutorsFER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧率

    % 为被攻击的执行体运行结果存在这
    TotalNoneAttackExecutorsBERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误码数
    % 例如：TotalExecutorsBERNum(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码数
    
    TotalNoneAttackExecutorsFERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误帧数
    % 例如：TotalExecutorsFERNum(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧数
    
    TotalNoneAttackExecutorsBER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误码率
    % 例如：TotalExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率
    
    TotalNoneAttackExecutorsFER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行完所有帧后的的误帧率
    % 例如：TotalExecutorsFER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧率
    
    EachFrameExecutorsBERs = zeros(8,total_frame_num); %存取每个执行体在各个帧下的误码率 最后3个是3余度dhr、4余度dhr、5余度dhr
    EachFrameNoneAttackExecutorsBERs = zeros(5,total_frame_num);% 无攻击的五个执行体的运行结果



    % 得到各执行体出错的帧数组
    SCBug = sort(randperm(total_frame_num , ceil(total_frame_num * SC_attack_rate)));
    BPBug = sort(randperm(total_frame_num , ceil(total_frame_num * BP_attack_rate)));
    SCANBug = sort(randperm(total_frame_num , ceil(total_frame_num * SCAN_attack_rate)));
    SCLBug = sort(randperm(total_frame_num , ceil(total_frame_num * SCL_attack_rate)));
    SSCBug = sort(randperm(total_frame_num , ceil(total_frame_num * SSC_attack_rate)));

    SCBug_index = 1;
    BPBug_index = 1;
    SCANBug_index = 1;
    SCLBug_index = 1;
    SSCBug_index =1;

    while IsContinue
        if current_total_frame_num>=total_frame_num
            break;
        else
%             ExecutorsBERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误码数
%             % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码数
%             
%             ExecutorsFERNum = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误帧数
%             % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧数
%             
%             ExecutorsBER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误码率
%             % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误码率
%             
%             ExecutorsFER = zeros(5,length(ebn0));%存取各个执行体在不同的信噪比条件下执行的误帧率
%             % 例如：ExecutorsBER(executors(1),ebn0(1)) 代表第一个执行体在第一个信噪比条件下执行后得到的误帧率

            executors = [1,2,3,4,5];%对比试验执行体集

            % 选取执行体构成执行体集合
%             executornums = [3,4,5]; %input('请输入需要选取的执行体数量（请请输入最小为1，最大为5的正整数）：');
%              
% %            for index = 1:length(executornums)
                [DHRExecutors,ExecutorConditions] = SwichExecutor(ExecutorConditions,redundancy);   
 

                for l = current_total_frame_num + 1 : current_total_frame_num + max_frame_num
                    
    %                 fprintf('\n Now running:%f  [%d of %d] \n\t Iteration-Counter: %53d',ebn0(j),j,length(ebn0),0);
    %                 tt=tic();%用于启动一个计时器以测量代码的执行时间。返回值tt是一个计时器对象，用于停止计时器并计算经过的时间 
                    
                    if attack_method == 1
                        % 每个帧选择被攻击的执行体数量
                        Attack_executor_num = randperm(redundancy,1);
                        % 选出被攻击的执行体集合
                        Attack_exexcutors = SwitchAttackExecutors(DHRExecutors, Attack_executor_num,redundancy);
%                         fprintf('\n 第%5d 帧被攻击执行体编号如下：\n',l);
%                         disp(Attack_exexcutors);
                    else
                        Attack_exexcutors = 0;
                    
                    end
    
                    %每次先生成信息序列，然后执行，使每个执行体执行相同的信息序列
                    u = randi(2,1,K)-1; %Bernoulli(0.5);  伯努利分布 生成一个长度为K的随机二进制序列（0或1），每个元素的取值都是等概率的
                    u_executors = zeros(5,K);
                    u_NoneAttackExecutors = zeros(5,K);
                    
                    % 生成无攻击时的信息
                    for i1 = 1:length(executors)
                        switch executors(i1)
                            case 1
%                                     fprintf('\n running polar_SC_decode now！');
                                crc_size = 0;
                                initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                            case 2
%                                     fprintf('\n running polar_BP_decode now！');
                                crc_size = 0;
                                initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                            case 3
%                                     fprintf('\n running polar_SCAN_decode now！');
                                crc_size = 0;
                                initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                            case 4
%                                     fprintf('\nrunning polar_SCL_decode now！');
                                initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                            case 5
%                                     fprintf('\n running polar_SSC_decode now！');
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
                        [TotalNoneAttackExecutorsBERNum(executors(i1),1),TotalNoneAttackExecutorsFERNum(executors(i1),1),u_NoneAttackExecutors(i1,:)] = ...
                            RunABlock(N,K,u,executors(i1),ebn0,TotalNoneAttackExecutorsBERNum(executors(i1),1),TotalNoneAttackExecutorsFERNum(executors(i1),1),...
                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);

                        TotalNoneAttackExecutorsBER(executors(i1),1) = TotalNoneAttackExecutorsBERNum(executors(i1),1)/(K*l);
                        TotalNoneAttackExecutorsFER(executors(i1),1) = TotalNoneAttackExecutorsFERNum(executors(i1),1)/l;
                        EachFrameNoneAttackExecutorsBERs(executors(i1),l) = TotalNoneAttackExecutorsBER(executors(i1),1);

                        

                         %计算误码数和误帧数
                        if l == current_total_frame_num + max_frame_num
                            fprintf(' \n无攻击时，异构执行体：%7s   在码率为：%7d ，信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误码数为：%7d ', ...
                                ExecutorNames(executors(i1)),K/N,ebn0,l-max_frame_num+1,l,TotalNoneAttackExecutorsBERNum(executors(i1),1));
                            fprintf(' \n无攻击时，异构执行体：%7s   在码率为：%7d ，信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误帧数为：%7d ', ...
                                ExecutorNames(executors(i1)),K/N,ebn0,l-max_frame_num+1,l,TotalNoneAttackExecutorsFERNum(executors(i1),1));
        %                             fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
                            fprintf('\n');
                        end
                    end
                      

    
                    %每个帧/块需要每个算法都执行一次（被攻击）
                    for ii = 1:length(executors)
                        
                        switch executors(ii)
                            case 1
                                % 判断是否被攻击
                                % 未被攻击
                                if isempty(SCBug)
                                    % 该执行体被选中
                                    if isempty(find(DHRExecutors == 1))
                                        SC_decoder_conditions(l) = 1;
                                    end
                                    crc_size = 0;
                                    initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                    %运行 得到误码数和误帧数
                                    [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                        RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                        bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                % 被攻击
                                else
                                    %逐帧分析是否被攻击
                                    if SCBug(SCBug_index) == l
                                        % 本次攻击力度
                                        SC_attack = roundn( 0.1 + (0.8 - 0.1)* rand(1) ,-1 );
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunAnAttackedBlock1(N,K,u,SC_attack,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
    %                                     ExecutorsBERNum(executors(ii),1) = ExecutorsBERNum(executors(ii),1) + ceil(SC_attack * K);
    %                                     ExecutorsFERNum(executors(ii),1) = ExecutorsFERNum(executors(ii),1) + 1;
                                       if SCBug_index < length(SCBug)
                                           SCBug_index = SCBug_index + 1;
                                           SC_decoder_conditions(l) = -1;
                                       end
                                    else
                                        % 该执行体被选中
                                        if isempty(find(DHRExecutors == 1))
                                            SC_decoder_conditions(l) = 1;
                                        end
                                        crc_size = 0;
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        %运行 得到误码数和误帧数
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                    end
                                end
                            case 2
                                % 判断是否被攻击
                                % 未被攻击
                                if isempty(BPBug)
                                    if isempty(find(DHRExecutors == 2))
                                        BP_decoder_conditions(l) = 1;
                                    end
                                    crc_size = 0;
                                    initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                    %运行 得到误码数和误帧数
                                    [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                        RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                        bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                % 被攻击
                                else
                                    if BPBug(BPBug_index) == l
                                        BP_attack = roundn(  0.1 + (0.8 - 0.1)* rand(1) ,-1 );
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunAnAttackedBlock1(N,K,u,BP_attack,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
    %                                     ExecutorsBERNum(executors(ii),1) = ExecutorsBERNum(executors(ii),1) + ceil(BP_attack * K);
    %                                     ExecutorsFERNum(executors(ii),1) = ExecutorsFERNum(executors(ii),1) + 1;
                                        if BPBug_index < length(BPBug)
                                            BPBug_index = BPBug_index + 1;
                                            BP_decoder_conditions(l) = -1;
                                        end
                                    else
                                        % 该执行体被选中
                                        if isempty(find(DHRExecutors == 2))
                                            BP_decoder_conditions(l) = 1;
                                        end
                                        crc_size = 0;
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        %运行 得到误码数和误帧数
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                    end
                                end
                            case 3
                                % 判断是否被攻击
                                % 未被攻击
                                if isempty(SCANBug)
                                    if isempty(find(DHRExecutors == 3))
                                        SCAN_decoder_conditions(l) = 1;
                                    end
                                    crc_size = 0; 
                                    initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                    %运行 得到误码数和误帧数
                                    [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                        RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                        bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                % 被攻击
                                else
                                    if SCANBug(SCANBug_index) == l
                                        SCAN_attack = roundn( 0.1 + (0.8 - 0.1)* rand(1) ,-1 );
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunAnAttackedBlock1(N,K,u,SCAN_attack,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
    %                                    ExecutorsBERNum(executors(ii),1) = ExecutorsBERNum(executors(ii),1) + ceil(SCAN_attack * K);
    %                                    ExecutorsFERNum(executors(ii),1) = ExecutorsFERNum(executors(ii),1) + 1;
                                       if SCANBug_index < length(SCANBug)
                                           SCANBug_index = SCANBug_index + 1;
                                           SCAN_decoder_conditions(l) = -1;
                                       end
                                    else
                                        % 该执行体被选中
                                        if isempty(find(DHRExecutors == 3))
                                            SCAN_decoder_conditions(l) = 1;
                                        end
                                        crc_size = 0; 
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        %运行 得到误码数和误帧数
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                    end
                                end
                            case 4
                                % 判断是否被攻击
                                % 未被攻击
                                if isempty(SCLBug)
                                    if isempty(find(DHRExecutors == 4))
                                        SCL_decoder_conditions(l) = 1;
                                    end
                                    initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                    %运行 得到误码数和误帧数
                                    [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                        RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                        bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                % 被攻击
                                else
                                    if SCLBug(SCLBug_index) == l
                                        SCL_attack = roundn( 0.1 + (0.8 - 0.1)* rand(1) ,-1 );
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunAnAttackedBlock1(N,K,u,SCL_attack,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
    %                                    ExecutorsBERNum(executors(ii),1) = ExecutorsBERNum(executors(ii),1) + ceil(SCL_attack * K);
    %                                    ExecutorsFERNum(executors(ii),1) = ExecutorsFERNum(executors(ii),1) + 1;
                                       if SCLBug_index < length(SCLBug)
                                           SCLBug_index = SCLBug_index + 1;
                                           SCL_decoder_conditions(l) = -1;
                                       end
                                    else
                                        % 该执行体被选中
                                        if isempty(find(DHRExecutors == 4))
                                            SCL_decoder_conditions(l) = 1;
                                        end
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        %运行 得到误码数和误帧数
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                    end
                                end
                            case 5
                                % 判断是否被攻击
                                % 未被攻击
                                if isempty(SSCBug)
                                    if isempty(find(DHRExecutors == 5))
                                        SSC_decoder_conditions(l) = 1;
                                    end
                                    crc_size = 0;
                                    initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                    [decoder_tree_initial, G_set, B_set] = intial_tree_G( );
                                    %运行 得到误码数和误帧数
                                    [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                        RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                        bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                % 被攻击
                                else
                                    if SSCBug(SSCBug_index) == l
                                        SSC_attack = roundn( 0.1 + (0.8 - 0.1)* rand(1) ,-1 );
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunAnAttackedBlock1(N,K,u,SSC_attack,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
    %                                    ExecutorsBERNum(executors(ii),1) = ExecutorsBERNum(executors(ii),1) + ceil(SSC_attack * K);
    %                                    ExecutorsFERNum(executors(ii),1) = ExecutorsFERNum(executors(ii),1) + 1;
                                       if SSCBug_index < length(SSCBug)
                                           SSCBug_index = SSCBug_index + 1;
                                           SSC_decoder_conditions(l) = -1;
                                       end
                                    else
                                        % 该执行体被选中
                                        if isempty(find(DHRExecutors == 5))
                                            SSC_decoder_conditions(l) = 1;
                                        end
                                        crc_size = 0;
                                        initPC(N,K,construction_method,design_snr_dB,sigma,crc_size);
    %                                     [decoder_tree_initial, G_set, B_set] = intial_tree_G( );
                                        %运行 得到误码数和误帧数
                                        [TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),u_executors(ii,:)] = ...
                                            RunABlock(N,K,u,executors(ii),ebn0,TotalExecutorsBERNum(executors(ii),1),TotalExecutorsFERNum(executors(ii),1),...
                                            bp_iter_num,scan_iter_num,scl_list_size,decoder_tree_initial, G_set, B_set);
                                     end
                                end
    
                            otherwise
                                fprintf('\n invalid input!!!');
                                bp_iter_num = 60;
                                scan_iter_num = 8;
                                scl_list_size = 4;
                                crc_size = 0;
                        end
                  
    
                        
                        
                        %计算误码数和误帧数
                        if l == current_total_frame_num + max_frame_num
                            fprintf(' \n随机攻击时，异构执行体：%7s   在码率为：%7d ，信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误码数为：%7d ',...
                                ExecutorNames(executors(ii)),K/N,ebn0,l-max_frame_num+1,l,TotalExecutorsBERNum(executors(ii),1));
                            fprintf(' \n随机攻击时，异构执行体：%7s   在码率为：%7d ，信噪比为：%7d  的环境下从第%7d 个块执行到第%7d 个块的总误帧数为：%7d ',...
                                ExecutorNames(executors(ii)),K/N,ebn0,l-max_frame_num+1,l,TotalExecutorsFERNum(executors(ii),1));
    %                             fprintf('\n\t Total time taken: %.2f sec (%d samples)',toc(tt),l);
                            fprintf('\n');
                        end
    
    
                        
                        %每次执行一个帧，计算每个执行体在不同信噪比下的BER、FER
    %                     if l == current_total_frame_num + max_frame_num
                            %求出每个执行体执行所有帧后的总BERNum、FERNum
    %                         TotalExecutorsBERNum(executors(ii),1) = TotalExecutorsBERNum(executors(ii),1);
    %                         TotalExecutorsFERNum(executors(ii),1) = TotalExecutorsFERNum(executors(ii),1);
    %                         ExecutorsBER(executors(ii),1) = TotalExecutorsBERNum(executors(ii),1)/(K*l);
    %                         ExecutorsFER(executors(ii),1) = TotalExecutorsFERNum(executors(ii),1)/l;
                            %得到每个执行体执行所有帧后的总BER、FER
                            TotalExecutorsBER(executors(ii),1) = TotalExecutorsBERNum(executors(ii),1)/(K*l);
                            TotalExecutorsFER(executors(ii),1) = TotalExecutorsFERNum(executors(ii),1)/l;
    %                     end
    
                        EachFrameExecutorsBERs(executors(ii),l) = TotalExecutorsBER(executors(ii),1);
                        
                    end
                    
     
                    
%                     switch redundancy
%                         case 3
%                             Three_u_executors = zeros(3,K);
%                             num3 = 1;
%                             for x = 1:length(executors)
%                                 if ismember(executors(x) ,DHRExecutors) == 1     
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1% 有攻击
%                                         Three_u_executors(num3,:) = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Three_u_executors(num3,:) = u_NoneAttackExecutors(executors(x),:);
%                                     end
% %                                 else
% %                                     Four_u_executors(executors(x),:) = u_NoneAttackExecutors(executors(x),:);
%                                     num3 = num3 + 1;
%                                 end
%                                
%                             end
%                             u_System = mode(Three_u_executors);
%                             %计算每一帧系统误码数
%                             nfails = sum(u_System ~= u);
%                             ThreeSystemBERNum = ThreeSystemBERNum + nfails;
%                             ThreeSystemFERNum = ThreeSystemFERNum + (nfails>0);
%                             ThreeSystemBER = ThreeSystemBERNum/(l*K);
%                             ThreeSystemFER = ThreeSystemFERNum/l;
%                             EachFrameExecutorsBERs(6,l) = ThreeSystemBER;
%                         case 4
%                             Four_u_executors = zeros(4,K);
%                             % 对每个执行体判断是否是执行体集合中的元素，是则进行下一步判断，不是则直接跳出
%                             % 下一步：判断执行体集合中的执行体是否被攻击，时则从u_executors中拿数据，不是则从u_NoneAttackExecutors中拿数据
%                             num4 = 1;
%                             for x = 1:length(executors)
%                                 if ismember(executors(x) ,DHRExecutors) == 1     
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1% 有攻击
%                                         Four_u_executors(num4,:) = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Four_u_executors(num4,:) = u_NoneAttackExecutors(executors(x),:);
%                                     end
% %                                 else
% %                                     Four_u_executors(executors(x),:) = u_NoneAttackExecutors(executors(x),:);
%                                     num4 = num4 + 1;
%                                 end
%                                 
%                             end
%                             u_System = mode(Four_u_executors);
%                             %计算每一帧系统误码数
%                             nfails = sum(u_System ~= u);
%                             FourSystemBERNum = FourSystemBERNum + nfails;
%                             FourSystemFERNum = FourSystemFERNum + (nfails>0);
%                             FourSystemBER = FourSystemBERNum/(l*K);
%                             FourSystemFER = FourSystemFERNum/l;
%                             EachFrameExecutorsBERs(7,l) = FourSystemBER;
%                         case 5
%                             Five_u_executors = zeros(5,K);
%                             num5 = 1;
%                             for x = 1:length(executors)
%                                 if ismember(executors(x) ,DHRExecutors) == 1     
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1% 有攻击
%                                         Five_u_executors(num5,:) = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Five_u_executors(num5,:) = u_NoneAttackExecutors(executors(x),:);
%                                     end
% %                                 else
% %                                     Four_u_executors(executors(x),:) = u_NoneAttackExecutors(executors(x),:);
%                                     num5 = num5 + 1;
%                                 end
%                                 
%                             end
%                             u_System = mode(Five_u_executors,1);
%                             %计算每一帧系统误码数
%                             nfails = sum(u_System ~= u);
%                             FiveSystemBERNum = FiveSystemBERNum + nfails;
%                             FiveSystemFERNum = FiveSystemFERNum + (nfails>0);
%                             FiveSystemBER = FiveSystemBERNum/(l*K);
%                             FiveSystemFER = FiveSystemFERNum/l;
%                             EachFrameExecutorsBERs(8,l) = FiveSystemBER;
%                     end
    
                    % 适用于无攻击时 对比不同冗余度的执行体
                    % 裁决 对比 u_executors和u
                    %3余度dhr
                    Three_u_executors = u_executors(DHRExecutors,:);
                    u_System = mode(Three_u_executors,1);
                    %计算每一帧系统误码数
                    nfails = sum(u_System ~= u);
                    ThreeSystemBERNum = ThreeSystemBERNum + nfails;
                    ThreeSystemFERNum = ThreeSystemFERNum + (nfails>0);
                    ThreeSystemBER = ThreeSystemBERNum/(l*K);
                    ThreeSystemFER = ThreeSystemFERNum/l;
                    EachFrameExecutorsBERs(6,l) = ThreeSystemBER;
    
                    %4余度dhr
                    %找到剩余执行体
                    temp_executors = executors( ~ismember(executors,DHRExecutors) );
                    temp_executors = temp_executors(2);
                    temp_executors = [DHRExecutors,temp_executors];
                    Four_u_executors = u_executors(temp_executors,:);
                    u_System = mode(Four_u_executors,1);
                    %计算每一帧系统误码数
                    nfails = sum(u_System ~= u);
                    FourSystemBERNum = FourSystemBERNum + nfails;
                    FourSystemFERNum = FourSystemFERNum + (nfails>0);
                    FourSystemBER = FourSystemBERNum/(l*K);
                    FourSystemFER = FourSystemFERNum/l;
                    EachFrameExecutorsBERs(7,l) = FourSystemBER;
    
                    %5余度dhr
                    u_System = mode(u_executors,1);
                    %计算每一帧系统误码数
                    nfails = sum(u_System ~= u);
                    FiveSystemBERNum = FiveSystemBERNum + nfails;
                    FiveSystemFERNum = FiveSystemFERNum + (nfails>0);
                    FiveSystemBER = FiveSystemBERNum/(l*K);
                    FiveSystemFER = FiveSystemFERNum/l;
                    EachFrameExecutorsBERs(8,l) = FiveSystemBER;
    
                    
    
                    
                
                end
                
    %             title('码率：',K/N);
    %             xlabel('已执行帧数')
    %             ylabel('误码率(BER)')
    %             legend('FER of DHR polar decoder','FER of SC polar decoder','FER of BP polar decoder' ,...
    %             'FER of SCAN polar decoder','FER of SCL polar decoder','FER of SSC polar decoder')
    %             hold on;
                
                %所有执行体执行完一轮（max_frame_num）后，比较被选出的执行体的BER和FER，并选最好的作为系统的结果,最差的下次不被选用
                [MINBERs,MINBERExecutors] = GetMinBER(length(DHRExecutors),TotalExecutorsBERNum,ebn0,DHRExecutors,ExecutorNames);
                [MINFERs,MINFERExecutors] = GetMinFER(length(DHRExecutors),TotalExecutorsFERNum,ebn0,DHRExecutors,ExecutorNames);
                [MAXBERs,MAXBERExecutors] = GetMaxBER(length(DHRExecutors),TotalExecutorsBERNum,ebn0,DHRExecutors,ExecutorNames);
                [MAXFERs,MAXFERExecutors] = GetMaxFER(length(DHRExecutors),TotalExecutorsFERNum,ebn0,DHRExecutors,ExecutorNames);
                fprintf(' \n异构执行体：%7s   不会在下一轮执行体选取中被选中 ',ExecutorNames(MAXBERExecutors));
                
                %ExecutorConditions全部清零，以便下次选取
                for i = 1:5
                    ExecutorConditions(i,1) = 0;
                end
    
                if MINBERExecutors ~= MAXBERExecutors %排除误码数为0时剔除一个异构体
                    ExecutorConditions(MAXBERExecutors,1) = 1;
                end
            
            %计算系统的误码总数和误帧总数
%             SystemBERNum = SystemBERNum + MINBERs;
%             SystemFERNum = SystemFERNum + MINFERs;
%             SystemBER = SystemBERNum/((current_total_frame_num+max_frame_num)*K);
%             SystemFER = SystemFERNum/(current_total_frame_num+max_frame_num);
    
    
                current_total_frame_num = current_total_frame_num + max_frame_num;
%             end
    
        end
       
    end





end

