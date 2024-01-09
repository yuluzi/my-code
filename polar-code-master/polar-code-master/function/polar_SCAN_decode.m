function [u_llr, c_llr] = polar_SCAN_decode(y_llr,iter_num)
% 这段 MATLAB 代码实现了极化码的 SCAN (Successive Cancellation Approximation Network) 译码算法，用于解码极化码的二元对称信道 (BSC) 上的传输信号。
% 
% 输入参数y_llr是接收到的二元对称信道上的对数似然比 (LLR) 信息。
% 输入参数iter_num是SCAN 迭代的次数。（数值结果越好但同时运行速度更慢）

% 算法首先根据极化码的码长 N 和信息长度 n，初始化左信息 L 和右信息 B 矩阵。
% 接下来进入主循环，对每个 phi，通过 updateLLRMap 函数更新左信息 L，
% 然后判断 phi 是否为偶数，如果不是，则通过 updateBitMap 函数更新右信息 B。
% 在循环结束后，输出最终的左信息和右信息，即 u_llr 和 c_llr。
% 
% 需要注意的是，在代码中有两行被注释掉的代码，这是 SCAN 算法中的一个变种，称为 SCAN-B，该变种在原有的 SCAN 算法的基础上加入了一些修正项。

    %初始化PCparams.L 和 PCparams.B
    global PCparams;
    N = PCparams.N;
    n = PCparams.n;
    
    plus_infinity = 1000;
    PCparams.L = zeros(N,n+1);%left message
    PCparams.B = zeros(N,n+1);%right message
    PCparams.L(:,n+1) = y_llr';%initial L
    PCparams.B(PCparams.FZlookup==0,1) = plus_infinity;%initial B

    
    %主循环
    for ii = 1:iter_num
        for phi = 0:N-1
            updateLLRMap(n,phi);

            if mod(phi,2)~=0
                %修订=============================
                 %PCparams.B(phi:phi+1,1) =  PCparams.B(phi:phi+1,1)+0.1*PCparams.L(phi:phi+1,1);
                 %================================
                updateBitMap(n,phi);          
            end
        end
        %PCparams.B =  PCparams.B+PCparams.L;
    end
       
    %输出最终的左信息u_llr和右信息c_llr
    u_llr = PCparams.L(:,1)+PCparams.B(:,1);
    
    c_llr = PCparams.B(:,n+1);
    
    
end
