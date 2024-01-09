function [u_llr, c_llr] = polar_SCAN_decode(y_llr,iter_num)
% ��� MATLAB ����ʵ���˼������ SCAN (Successive Cancellation Approximation Network) �����㷨�����ڽ��뼫����Ķ�Ԫ�Գ��ŵ� (BSC) �ϵĴ����źš�
% 
% �������y_llr�ǽ��յ��Ķ�Ԫ�Գ��ŵ��ϵĶ�����Ȼ�� (LLR) ��Ϣ��
% �������iter_num��SCAN �����Ĵ���������ֵ���Խ�õ�ͬʱ�����ٶȸ�����

% �㷨���ȸ��ݼ�������볤 N ����Ϣ���� n����ʼ������Ϣ L ������Ϣ B ����
% ������������ѭ������ÿ�� phi��ͨ�� updateLLRMap ������������Ϣ L��
% Ȼ���ж� phi �Ƿ�Ϊż����������ǣ���ͨ�� updateBitMap ������������Ϣ B��
% ��ѭ��������������յ�����Ϣ������Ϣ���� u_llr �� c_llr��
% 
% ��Ҫע����ǣ��ڴ����������б�ע�͵��Ĵ��룬���� SCAN �㷨�е�һ�����֣���Ϊ SCAN-B���ñ�����ԭ�е� SCAN �㷨�Ļ����ϼ�����һЩ�����

    %��ʼ��PCparams.L �� PCparams.B
    global PCparams;
    N = PCparams.N;
    n = PCparams.n;
    
    plus_infinity = 1000;
    PCparams.L = zeros(N,n+1);%left message
    PCparams.B = zeros(N,n+1);%right message
    PCparams.L(:,n+1) = y_llr';%initial L
    PCparams.B(PCparams.FZlookup==0,1) = plus_infinity;%initial B

    
    %��ѭ��
    for ii = 1:iter_num
        for phi = 0:N-1
            updateLLRMap(n,phi);

            if mod(phi,2)~=0
                %�޶�=============================
                 %PCparams.B(phi:phi+1,1) =  PCparams.B(phi:phi+1,1)+0.1*PCparams.L(phi:phi+1,1);
                 %================================
                updateBitMap(n,phi);          
            end
        end
        %PCparams.B =  PCparams.B+PCparams.L;
    end
       
    %������յ�����Ϣu_llr������Ϣc_llr
    u_llr = PCparams.L(:,1)+PCparams.B(:,1);
    
    c_llr = PCparams.B(:,n+1);
    
    
end
