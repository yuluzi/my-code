function u_llr = polar_SC_decode(initial_llr)
% 将初始化 LLR（对数似然比）作为输入，并返回解码器的输出 u_llr。
% 
% 主函数中，首先从全局变量 PCparams 中获取码长 N 和信息位长度 n 的信息，
% 并初始化两个大小为 N x (n+1) 的矩阵 L 和 B。
% 然后，将输入的初始化 LLR 存储在 L 矩阵的最后一列中。
% 之后，从 0 到 N-1 的枚举变量 phi 循环，分别调用 updateL 和 updateB 函数对矩阵 L 和 B 进行更新。
% 在更新完成后，将 L 矩阵的第一列作为输出 u_llr 并返回。
% 
% updateL 函数实现了 SC 译码算法的主要部分。
% 在该函数中，首先检查是否达到了最后一层，如果达到最后一层，则直接返回。
% 否则，计算中间变量 psi 并根据 phi 的奇偶性递归调用 updateL 函数。
% 在递归完成后，针对所有的 omega（0 到 2^(n-lambda)-1），根据 phi 的奇偶性更新 L 矩阵的值。
% 如果 phi 是偶数，则根据算法中的定义计算 fFunction 并存储到 L 矩阵中；
% 否则，根据 B 矩阵中的值计算并更新 L 矩阵中的值。
% 
% updateB 函数用于更新 B 矩阵中的值。
% 在该函数中，首先计算中间变量 psi 并检查是否到达最后一层。
% 如果未到达最后一层，则递归调用 updateB 函数。
% 在递归完成后，对于所有 omega（0 到 2^(n-lambda)-1），根据算法中的定义更新 B 矩阵中的值。


    global PCparams;
    N = PCparams.N;
    n = PCparams.n;
    PCparams.L = zeros(N,n+1);
    PCparams.B = zeros(N,n+1);
    PCparams.L(:,n+1) = initial_llr';
    
    for phi = 0:N-1
        updateL(n,phi);
        if PCparams.FZlookup(phi+1) == 0
            PCparams.B(phi+1,1) =  0;
        else
            if PCparams.L(phi+1,1)<0
                PCparams.B(phi+1,1) = 1;
            else
                PCparams.B(phi+1,1) = 0;
            end
        end
        if mod(phi,2)==1
            updateB(n,phi);
        end
    end
    
    u_llr = PCparams.L(:,1);
    
end


% function updateL(lambda,phi)
%     global PCparams;
%     n = PCparams.n;
%     if lambda == 0
%         return;
%     end
%     psi = floor(phi/2);
%     if mod(phi,2)==0
%         updateL(lambda-1,psi);
%     end
%     for omega=0:2^(n-lambda)-1
%         if mod(phi,2)==0
%             %do sth
%             PCparams.L(phi+omega*2^lambda+1,n+1-lambda) = fFunction(PCparams.L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),PCparams.L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda));
%         else
%             %do sth
%             if PCparams.B(phi-1+omega*2^(lambda)+1,n+1-lambda) == 0
%                 PCparams.L(phi+omega*2^(lambda)+1,n+1-lambda) = PCparams.L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+PCparams.L(psi+(2*omega)*2^(lambda-1)+1,n+2-lambda);
%             else
%                 PCparams.L(phi+omega*2^(lambda)+1,n+1-lambda) = PCparams.L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)-PCparams.L(psi+(2*omega)*2^(lambda-1)+1,n+2-lambda);
%             end
%         end
%     end
% end

% function updateB(lambda,phi)
%     global PCparams;
%     n = PCparams.n;
%     
%     psi = floor(phi/2);
%     if mod(phi,2)~=0
%         for omega = 0:2^(n-lambda)-1
%             PCparams.B(psi+2*omega*2^(lambda-1)+1,n+2-lambda) = xor(PCparams.B(phi-1+omega*2^(lambda)+1,n+1-lambda),PCparams.B(phi+omega*2^(lambda)+1,n+1-lambda));
%             PCparams.B(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda) = PCparams.B(phi+omega*2^(lambda)+1,n+1-lambda);
%         end
%         if mod(psi,2)~=0
%             updateB(lambda-1,psi);
%         end
%     end
% end