%使用巴氏参数上界方法构造polar code
%polar encoder and SC, BP , SCAN decoder for awgn channel
%panzhipeng
function construct_polar_code_Ba(N,design_snr_dB)

    %计算极化码的阶数n，并检查输入参数N是否是2的幂次方。如果不是，则输出一条错误信息并结束函数。
    n = ceil(log2(N)); 
    NN = 2^n;
    if(NN~=N)
        fprintf('The num N must be the power of 2!');
        return;
    end

    %创建一个文件，文件名包含码长和设计信噪比
    file_name = sprintf('PolarCode_block_length_%d_designSNR_%.2fdB_method_BhattaBound.txt',N,design_snr_dB);
    fid = fopen(file_name,'w+');

    %计算比特反转置换表，将所有索引值按比特反转的顺序排列。
    bitreversedindices = zeros(1,N);
    for index = 1 : N
        bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,n)));
    end
    
    %计算巴氏参数上界。首先，初始化z(1)为设计信噪比的巴氏参数上界。
    %然后，按照从低到高的层次，使用递归的方式计算巴氏参数上界。
    z = zeros(1,N);
    design_snr_num = 10^(design_snr_dB/10);
    z(1) = log(exp(-design_snr_num));
    %z(1) = -1;%designSNR = 0dB, z(1) = exp(-S) in logdomain z(1) = -designSNR_num = -1
    for lev = 1:n
        B = 2^lev;
        for j = 1:B/2
            T = z(j);
            z(j) = log(2)+T + log1p(-exp(T-log(2)));
            z(B/2+j) = 2*T;
        end
    end
    
    %将z值按比特反转的顺序排列，并对排列后的索引进行排序。
    %将排序后的索引写入文件中，结束函数。这些索引用于后续的编码和解码过程
    z = z(bitreversedindices+1);
    [~,indices] = sort(z,'ascend');
    for ii = 1:length(z)
        fprintf(fid,'%d\r\n',indices(ii));
    end
    fclose(fid);
end
