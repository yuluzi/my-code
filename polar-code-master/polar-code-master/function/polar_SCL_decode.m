function [u_llr] = polar_SCL_decode(llr, list_size)
% 该MATLAB函数实现了极化分组码的逐步递归译码（SCL译码）算法。
% 函数的输入包括信息位的对数似然比（LLR）和列表大小（list_size，数值结果越好但同时运行速度更慢）。
% 函数的输出是解码的信息位的对数似然比（u_llr）。
% 
% 该函数的算法流程如下：
% 初始化数据结构（initializeDataStructures），包括数组、链表、指针等。其中链表用于存储路径，指针用于快速访问链表中的元素。
% 分配初始路径（assignInitialPath），将所有路径初始化为0。函数返回初始路径的索引。
% 为初始路径分配空间，用llr更新llr_scl数组。 
% 对于每个极化位置phi=0,1,…,N-1进行如下操作： 
% a. 递归计算P_scl数组的值，其中P_scl用于存储每个路径的概率值。 
% b. 如果phi对应的位置是一个冻结位（FZlookup(phi+1)=0），则只保留路径上的冻结位（continuePaths_FrozenBit）；
%    否则，保留路径上的冻结位和非冻结位（continuePaths_UnfrozenBit）。
% c. 如果phi是奇数，则递归更新C_scl数组的值，其中C_scl用于存储每个路径的估计值。 
% 找到最可能的路径（findMostProbablePath），返回最可能的路径的索引。 
% 根据最可能的路径，返回信息位的对数似然比（u_llr）。 
% 值得注意的是，此函数的输入是极化分组码的对数似然比（LLR），而不是二进制数据。
% 因此，在第5步中找到最可能的路径时，需要将信息位从硬判决转换为软判决。



    global PCparams;
    PCparams.list_size = list_size;
    
    initializeDataStructures();
    l_index = assignInitialPath();
    s_index = getArrayPointer_P(0, l_index);
    PCparams.llr_scl( (get_i_scl(0, 0, s_index) + 1) : (get_i_scl(0, PCparams.N - 1, s_index) + 1) ) = llr;
    
    for phi = 0 : PCparams.N - 1
        
        recursivelyCalcP_scl(PCparams.n, phi);
        
        if PCparams.FZlookup(phi + 1) == 0
            continuePaths_FrozenBit(phi);
        else
            continuePaths_UnfrozenBit(phi);
        end
        
        if mod(phi, 2) == 1
            recursivelyUpdateC_scl(PCparams.n, phi);
        end
        
    end
    
    l_index = findMostProbablePath(1);
    c_m = getArrayPointer_C(PCparams.n, l_index);
    info = PCparams.i_scl(c_m+1,:);
    %u = info(PCparams.FZlookup == -1);
    %in order to make the function interface with the main
    %function, we let the hard decision of the info change to the
    %soft decision. note that this soft decision is not really the
    %true soft decision.
    u_llr = (1-2*info);
end