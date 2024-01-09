function [u_llr,c_llr] = polar_BP_decode(initialLLRs,BP_ITER_NUM)
%��Matlab����ʵ���˼������BP�����㷨�����У�
% �������initialLLRs�Ǽ�����ĳ�ʼ������Ȼ�ȣ�LLR����
% �������BP_ITER_NUMָ����BP�����㷨�ĵ�������������ֵ���Խ�õ�ͬʱ�����ٶȸ�����
% �������u_llr�ǽ������Ϣ���ض�Ӧ��LLR��
% �������c_llr�ǽ����У����ض�Ӧ��LLR��

% �ú�������Ҫ����������Ƕ�׵�forѭ�������е�һ��ѭ���Ա���lambda���е�����
% ���ڶ���ѭ���Ա���phi���е�����������ѭ�����������ʵ����BP�����㷨�е���Ϣ���ݹ��̡�
% ������˵������ÿ��lambda��phi�����㷨���������򴫵ݺ����򴫵ݵĸ�����Ϣ����������洢����������L��R�С�
% Ȼ�󣬸��㷨����lambda��phi��ֵ�����ظ�����BP��Ϣ���ݣ�ֱ���ﵽָ���ĵ���������

% ��ÿ�������Ľ���ʱ�����㷨������Ϣ���غ�У����ص�LLR����������Ϊ������ء�
% ��Ҫע����ǣ��ú�����ʹ����һ����ΪfFunction���Ӻ��������������Ϣ��
% ����Ӻ�����ʵ��ͨ����ͨ�����ʵ�ֵģ������㷨��ѡ�񣬿���ʹ�ò�ͬ�Ĳ����ԡ�

    global PCparams;
    n = PCparams.n;
    N = PCparams.N;
   
    L = zeros(N,n+1);
    R = zeros(N,n+1);
    inf_num = 1000;
    L(:,n+1) = initialLLRs';
    R(PCparams.FZlookup==0,1) = inf_num;
    
    for iter = 1:BP_ITER_NUM
        
        for lambda = 1:n
            for phi = 0:2^lambda-1
                psi = floor(phi/2);
				for omega = 0:2^(n-lambda)-1
					if mod(phi,2)==0	
						L(phi+omega*2^lambda+1,n+1-lambda)= fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+R(phi+1+omega*2^lambda+1,n+1-lambda));
					else
						L(phi+omega*2^lambda+1,n+1-lambda)= L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),R(phi-1+omega*2^lambda+1,n+1-lambda));
					end
				end
            end
        end
        
        for lambda = n:-1:1
            for phi = 0:2^lambda-1
				psi = floor(phi/2);
                if mod(phi,2)~=0	
                    for omega = 0:2^(n-lambda)-1
						R(psi+2*omega*2^(lambda-1)+1,n+2-lambda)= fFunction(L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+R(phi+omega*2^lambda+1,n+1-lambda), R(phi-1+omega*2^lambda+1,n+1-lambda));
						R(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)= R(phi+omega*2^lambda+1,n+1-lambda)+fFunction(L(psi+(2*omega)*2^(lambda-1)+1,n+2-lambda),R(phi-1+omega*2^lambda+1,n+1-lambda));
					end
			    end
            end
        end
        
    end
    
    %计算�?��输出llr
    u_llr = L(:,1)+R(:,1);
    c_llr = R(:,n+1);
end