%ʹ�ð��ϲ����Ͻ緽������polar code
%polar encoder and SC, BP , SCAN decoder for awgn channel
%panzhipeng
function construct_polar_code_Ba(N,design_snr_dB)

    %���㼫����Ľ���n��������������N�Ƿ���2���ݴη���������ǣ������һ��������Ϣ������������
    n = ceil(log2(N)); 
    NN = 2^n;
    if(NN~=N)
        fprintf('The num N must be the power of 2!');
        return;
    end

    %����һ���ļ����ļ��������볤����������
    file_name = sprintf('PolarCode_block_length_%d_designSNR_%.2fdB_method_BhattaBound.txt',N,design_snr_dB);
    fid = fopen(file_name,'w+');

    %������ط�ת�û�������������ֵ�����ط�ת��˳�����С�
    bitreversedindices = zeros(1,N);
    for index = 1 : N
        bitreversedindices(index) = bin2dec(wrev(dec2bin(index-1,n)));
    end
    
    %������ϲ����Ͻ硣���ȣ���ʼ��z(1)Ϊ�������ȵİ��ϲ����Ͻ硣
    %Ȼ�󣬰��մӵ͵��ߵĲ�Σ�ʹ�õݹ�ķ�ʽ������ϲ����Ͻ硣
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
    
    %��zֵ�����ط�ת��˳�����У��������к��������������
    %������������д���ļ��У�������������Щ�������ں����ı���ͽ������
    z = z(bitreversedindices+1);
    [~,indices] = sort(z,'ascend');
    for ii = 1:length(z)
        fprintf(fid,'%d\r\n',indices(ii));
    end
    fclose(fid);
end
