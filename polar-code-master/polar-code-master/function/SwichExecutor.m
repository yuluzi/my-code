function [executors,ExecutorConditions] = SwichExecutor(ExecutorConditions,n)
% 输入：
% ExecutorConditions:执行体状态数组 
% n:产生几个不同的异构执行体
% 输出：
% 输出产生的异构执行体集合executors
executors = zeros(1,n);
k = 1;
if n>5
    fprintf('\n 选取的异构执行体数量最大为5，请重新输入：\n');
    return
else
    %   随机产生一个当前不再执行状态的异构执行体编号
    fprintf('\n 选取的异构执行体编号如下：\n');
    while true
        if n > 0 %没选到n个执行体
            executor =  randperm(5,1);%产生1个1-5之间的随机数
            if ExecutorConditions(executor,1) == 0
                executors(1,k) = executor;
                k = k + 1;
                n = n - 1;
                ExecutorConditions(executor,1) = 1;
                disp(executor);
            else
                continue;
            end
        else
            break;
        end

    end
end



