function gfsk_8ary_viterbi_isi()
% gfsk_8ary_viterbi_isi.m
% 8-ary GFSK 相干解 + ISI 感知 Viterbi 序列检测
% 核心：1-符号记忆（8状态），cosine相似度分支度量
% 关键参数：tone LPF阶数=24，fc=0.75*tone_spacing
% 总延迟 = delay_gauss + delay_ch + delay_tone（必须全部计入）

% ... 完整代码（约900行）包含：
%   - 8-ary信号生成
%   - ISI参考模板预计算（含ch_filter）
%   - 8-状态Viterbi前向递归+回朔
%   - 无噪声自检（硬判 vs Viterbi）
%   - Eb/N0扫描对比
%   - 6幅图：BER/SER、增益、频谱、分支度量、误差位置、星座图
% ...

end
