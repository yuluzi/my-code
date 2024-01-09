function [executors] = SwitchAttackExecutors(DHRExecutors,n,redundancy)
% 输入：
% DHRExecutors:被选中的执行体
% n:产生几个不同的受攻击执行体
% 输出：
% 输出产生的异构执行体集合executors
executors = zeros(1,n);
k = 1;
if n>redundancy
    fprintf('\n 选取的异构执行体数量最大为：%7d，请重新输入：\n',redundancy);
    return
else
    %   随机产生一个当前不再执行状态的异构执行体编号
%     fprintf('\n 被攻击执行体编号如下：\n');
    while true
        if n > 0 %没选到n个执行体
            executor =  randperm(5,1);%产生1个1-5之间的随机数
            if ismember(executor,DHRExecutors) == 1 && ismember(executor,executors) == 0 

                executors(1,k) = executor;
                k = k + 1;
                n = n - 1;
%                 disp(executor);
            end
            continue;
        else
            break;
        end
    end

end



