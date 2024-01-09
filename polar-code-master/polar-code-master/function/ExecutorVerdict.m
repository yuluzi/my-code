function [outputArg1,outputArg2] = ExecutorVerdict(inputArg1,inputArg2)
%EXECUTORVERDICT 此处显示有关此函数的摘要
%   此处显示详细说明
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

 switch redundancy
%                         case 3
%                             for x = 1:length(DHRExecutors)
%                                 if attack_method == 1
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1 % 有攻击
%                                         Three_u_executors = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Three_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                     end
%                                 else
%                                     Three_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                 end
%                             end
%                             u_System = mode(Three_u_executors,1);
%                             %计算每一帧系统误码数
%                             nfails = sum(u_System ~= u);
%                             ThreeSystemBERNum = ThreeSystemBERNum + nfails;
%                             ThreeSystemFERNum = ThreeSystemFERNum + (nfails>0);
%                             ThreeSystemBER = ThreeSystemBERNum/(l*K);
%                             ThreeSystemFER = ThreeSystemFERNum/l;
%                             EachFrameExecutorsBERs(6,l) = ThreeSystemBER;
%                         case 4
%                             for x = 1:length(DHRExecutors)
%                                 if attack_method == 1
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1% 有攻击
%                                         Four_u_executors = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Four_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                     end
%                                 else
%                                     Four_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                 end
%                             end
%                             u_System = mode(Four_u_executors,1);
%                             %计算每一帧系统误码数
%                             nfails = sum(u_System ~= u);
%                             FourSystemBERNum = FourSystemBERNum + nfails;
%                             FourSystemFERNum = FourSystemFERNum + (nfails>0);
%                             FourSystemBER = FourSystemBERNum/(l*K);
%                             FourSystemFER = FourSystemFERNum/l;
%                             EachFrameExecutorsBERs(7,l) = FourSystemBER;
%                         case 5
%                             for x = 1:length(DHRExecutors)
%                                 if attack_method == 1
%                                     if ismember(executors(x) ,Attack_exexcutors) == 1% 有攻击
%                                         Five_u_executors = u_executors(executors(x),:);
%                                     else
%                                     % 无攻击
%                                         Five_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                     end
%                                 else
%                                     Five_u_executors = u_NoneAttackExecutors(executors(x),:);
%                                 end
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
    
                    %裁决 对比 u_executors和u
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