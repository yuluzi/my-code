function [MAXBERs,MAXBERExecutors] = GetMaxBER(n,ExecutorsBER,ebn0,executors,ExecutorNames)
% 本函数用于求出每个执行体集合在某个特定的信噪比环境下的执行后的误码率最大的执行体及其误码率
% 
% 输入：
% n：组成该执行体集合的个执行体数量 
% ExecutorsBER：各个执行体在各信噪比环境下执行得到的误码数累加和数组 
% ebn0：信噪比数组 
% executors：当前执行体集合 
% ExecutorNames：执行体名称数组 
% 
% 输出：
% MINBERs：该执行体集合在各信噪比环境下执行后误码率最大的执行体的误码率
% MINBERExecutors：该执行体集合在各信噪比环境下执行后误码率最大的执行体
% 


% 初始化每个执行体在每个信噪比环境下的最大BER或FER，MAXNums(1,1)代表第一个执行体在第一个信噪比下的最大BER或FER
MAXBERs = zeros(1,length(ebn0));

MAXBERExecutors = zeros(1,length(ebn0));
WorstExecutor = executors(1);

% 每个信噪比条件下对比每个执行体的BER
for j = 1:length(ebn0)
    max = -1;
    for i = 1:n
        if max<ExecutorsBER(executors(i),j)
            max = ExecutorsBER(executors(i),j);
            WorstExecutor = executors(i);
        end
    MAXBERs(1,j) = max;
    MAXBERExecutors(1,j) = WorstExecutor;
    end
    
end


for j = 1:length(ebn0)
    fprintf(' \n异构执行体：%7s   在信噪比为：%7d  的环境下总误码数最高，为：%7d ',ExecutorNames(MAXBERExecutors(1,j)),ebn0(j),MAXBERs(1,j));
end