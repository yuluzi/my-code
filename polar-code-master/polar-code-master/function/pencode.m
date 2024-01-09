function y=pencode(u, F_kron_n) 
% 输入：   
% u：长度为K的二进制向量，表示待编码的消息。
% F_kron_n：表示克罗内克矩阵的一个矩阵。
% 输出：
% y：长度为N的二进制向量，表示经过极化编码后的码字。

% 函数实现：
% 首先，代码读取全局变量PCparams中的FZlookup、crc_size和bitreversedindices参数。
% 然后，如果crc_size不为零，则执行循环冗余检查（CRC）编码，并将CRC编码后的结果追加到原始消息向量u的末尾。
% 接下来，用消息和CRC编码的结果替换FZlookup中值为-1的位置。
% 然后，将FZlookup中的元素重新排序，以适应PCparams.bitreversedindices中指定的位翻转顺序。
% 最后，将结果乘以克罗内克矩阵F_kron_n得到极化编码后的码字y，并返回结果。
global PCparams;

% FZlookup变量是一个长度为2^K的向量，表示在极化编码中所有长度为K的信息比特的所有可能输入位的极化转移序列。
% 其中，FZlookup中的值为-1表示该位置为校验位，需要通过CRC校验。
x = PCparams.FZlookup; 

if PCparams.crc_size ~=0
    
	crc = mod(PCparams.crc_matrix*u', 2)';
else
	crc = [];
end
u = [u crc];
% fprintf('x (x == -1)有:%3d位  ，u有:%3d位  \n',length(x (x == -1)),length(u));

x (x == -1) = u; % -1's will get replaced by message bits below
x = x(PCparams.bitreversedindices+1);
y = mod(x*F_kron_n,2);

end