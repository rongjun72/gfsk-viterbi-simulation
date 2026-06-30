import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================================
# 1. 参数配置
# ============================================================
MODULATION = 'BPSK'
M = 2
SNR_DB_RANGE = np.arange(0, 16, 2)
N_SYMBOLS = 20000
SEED = 42

# 3径 ISI 信道
CHANNEL_TAPS = np.array([0.8, 0.5, 0.3], dtype=complex)
L = len(CHANNEL_TAPS)

np.random.seed(SEED)

# ============================================================
# 2. BPSK 调制/解调
# ============================================================
def bpsk_modulate(bits):
    return 1 - 2 * bits.astype(float)

def bpsk_demodulate(symbols):
    return (np.real(symbols) < 0).astype(int)

# ============================================================
# 3. 信道仿真
# ============================================================
def transmit(symbols, h, snr_db):
    N = len(symbols)
    r_clean = np.convolve(symbols, h, mode='full')[:N]
    signal_power = np.mean(np.abs(r_clean) ** 2)
    snr_linear = 10 ** (snr_db / 10)
    noise_power = signal_power / snr_linear
    noise = np.sqrt(noise_power / 2) * (np.random.randn(N) + 1j * np.random.randn(N))
    return r_clean + noise

# ============================================================
# 4. Viterbi MLSE 均衡器（修正版）
# ============================================================
class ViterbiMLSE:
    """
    状态定义：state = [b_{k-1}, b_{k-2}, ..., b_{k-L+1}] 的位编码
      bits[0] = b_{k-1} (最新) → (state >> (L-2)) & 1
      bits[L-2] = b_{k-L+1} (最旧) → state & 1
    
    输入 s_k (bit=input)，新状态 = [b_k, b_{k-1}, ..., b_{k-L+2}]
      next_state = (input << (L-2)) | (state >> 1)
    
    期望输出 = h[0]*s_k + h[1]*s_{k-1} + ... + h[L-1]*s_{k-L+1}
    """
    def __init__(self, h):
        self.h = np.array(h, dtype=complex)
        self.L = len(self.h)
        self.M = 2
        self.n_states = self.M ** (self.L - 1)
        self.sym_map = np.array([1, -1], dtype=complex)
        self._build_trellis()
    
    def _build_trellis(self):
        self.trans = {}
        for s in range(self.n_states):
            bits = [(s >> (self.L - 2 - i)) & 1 for i in range(self.L - 1)]
            syms = self.sym_map[bits]
            for inp in [0, 1]:
                inp_sym = self.sym_map[inp]
                full = np.concatenate([[inp_sym], syms])
                expected = np.dot(self.h, full)
                next_s = (inp << (self.L - 2)) | (s >> 1)
                self.trans[(s, inp)] = (next_s, expected)
    
    def detect(self, r):
        N = len(r)
        INF = 1e18
        pm = np.full(self.n_states, INF)
        pm[0] = 0.0
        survivor = np.zeros((N, self.n_states), dtype=int) - 1
        
        for t in range(N):
            npm = np.full(self.n_states, INF)
            for s in range(self.n_states):
                if pm[s] >= INF:
                    continue
                for inp in [0, 1]:
                    ns, exp = self.trans[(s, inp)]
                    bm = np.abs(r[t] - exp) ** 2
                    cand = pm[s] + bm
                    if cand < npm[ns]:
                        npm[ns] = cand
                        survivor[t, ns] = s
            pm = npm
        
        state = int(np.argmin(pm))
        bits = np.zeros(N, dtype=int)
        for t in range(N - 1, -1, -1):
            ps = survivor[t, state]
            if ps < 0:
                break
            inp = (state >> (self.L - 2)) & 1
            bits[t] = inp
            state = ps
        return bits

# ============================================================
# 5. MMSE 线性均衡器
# ============================================================
class MMSEEqualizer:
    def __init__(self, h, snr_db, eq_len=11):
        self.h = np.array(h, dtype=complex)
        self.L = len(self.h)
        self.eq_len = eq_len
        self.snr = 10 ** (snr_db / 10)
        self._design()
    
    def _design(self):
        L, N = self.L, self.eq_len
        H = np.zeros((N + L - 1, N), dtype=complex)
        for i in range(N):
            H[i:i+L, i] = self.h
        R = H.T.conj() @ H + (1.0 / self.snr) * np.eye(N)
        p = H[L - 1, :].conj()
        self.w = np.linalg.solve(R, p)
    
    def equalize(self, r):
        return np.convolve(r, self.w, mode='valid')

# ============================================================
# 6. 运行仿真
# ============================================================
def main(ctx):
    out = []
    out.append("=" * 60)
    out.append("  最大似然序列估计 (MLSE) — Viterbi 算法演示")
    out.append("=" * 60)
    out.append(f"  调制方式: {MODULATION}")
    out.append(f"  信道抽头: {CHANNEL_TAPS}")
    out.append(f"  信道长度 L: {L}")
    out.append(f"  每 SNR 发送符号: {N_SYMBOLS}")
    out.append(f"  SNR 范围: {SNR_DB_RANGE[0]} ~ {SNR_DB_RANGE[-1]} dB")
    
    viterbi = ViterbiMLSE(CHANNEL_TAPS)
    out.append(f"\n  Viterbi 状态数: {viterbi.n_states} = 2^{L-1}")
    
    # Trellis 转移表
    out.append("\n" + "=" * 60)
    out.append("  Trellis 状态转移表")
    out.append("=" * 60)
    out.append(f"  {'状态':>6} {'符号历史':>12} | {'输入':>6} | {'下一状态':>8} | {'期望输出':>10}")
    out.append("  " + "-" * 50)
    for s in range(viterbi.n_states):
        bits = [(s >> (L - 2 - i)) & 1 for i in range(L - 1)]
        hist = "".join(["+" if b == 0 else "-" for b in bits])
        for inp in [0, 1]:
            ns, exp = viterbi.trans[(s, inp)]
            out.append(f"  {s:>6} {hist:>12} | {'+1' if inp == 0 else '-1':>6} | {ns:>8} | {exp.real:>10.3f}")
    
    ber_mls = []
    ber_mmse = []
    
    for snr_db in SNR_DB_RANGE:
        tx_bits = np.random.randint(0, 2, N_SYMBOLS)
        tx_symbols = bpsk_modulate(tx_bits)
        rx = transmit(tx_symbols, CHANNEL_TAPS, snr_db)
        
        # ---- MLSE ----
        prefix = np.zeros(L - 1, dtype=complex)
        ext_rx = np.concatenate([prefix, rx])
        det_mls = viterbi.detect(ext_rx)
        det_mls = det_mls[L - 1 : L - 1 + N_SYMBOLS]
        err_mls = int(np.sum(det_mls != tx_bits))
        ber_mls.append(float(err_mls / N_SYMBOLS))
        
        # ---- MMSE ----
        eq = MMSEEqualizer(CHANNEL_TAPS, snr_db, eq_len=11)
        rx_pad = np.concatenate([rx, np.zeros(eq.eq_len - 1, dtype=complex)])
        eq_out = eq.equalize(rx_pad)[:N_SYMBOLS]
        det_mmse = bpsk_demodulate(eq_out)
        err_mmse = int(np.sum(det_mmse != tx_bits))
        ber_mmse.append(float(err_mmse / N_SYMBOLS))
        
        out.append(f"\n  SNR = {snr_db:2d} dB:")
        out.append(f"    MLSE BER = {err_mls}/{N_SYMBOLS} = {err_mls/N_SYMBOLS:.2e}")
        out.append(f"    MMSE BER = {err_mmse}/{N_SYMBOLS} = {err_mmse/N_SYMBOLS:.2e}")
    
    # ============================================================
    # 7. 绘图
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    ax = axes[0]
    ax.semilogy(SNR_DB_RANGE, ber_mls, 'o-', linewidth=2.5, markersize=8,
                label='MLSE (Viterbi)', color='#1f77b4')
    ax.semilogy(SNR_DB_RANGE, ber_mmse, 's--', linewidth=2.5, markersize=8,
                label='MMSE Linear Equalizer', color='#ff7f0e')
    ax.set_xlabel('SNR (dB)', fontsize=12)
    ax.set_ylabel('Bit Error Rate (BER)', fontsize=12)
    ax.set_title('BER Performance: MLSE vs MMSE', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='lower left')
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.set_ylim([1e-4, 1])
    
    ax = axes[1]
    taps = np.arange(len(CHANNEL_TAPS))
    ax.stem(taps, np.abs(CHANNEL_TAPS), basefmt=' ', linefmt='C0-', markerfmt='C0o')
    ax.set_xlabel('Tap Index', fontsize=12)
    ax.set_ylabel('Magnitude', fontsize=12)
    ax.set_title('Channel Impulse Response (ISI Channel)', fontsize=14, fontweight='bold')
    ax.set_xticks(taps)
    ax.grid(True, alpha=0.3)
    for i, val in enumerate(np.abs(CHANNEL_TAPS)):
        ax.annotate(f'{val:.2f}', xy=(i, val), xytext=(5, 5),
                    textcoords='offset points', fontsize=10)
    
    plt.tight_layout()
    save_path = os.path.join(ctx["runDir"], 'mlse_vs_mmse_ber.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    import shutil
    ws_path = os.path.join("C:\\Users\\Administrator\\Documents\\Kimi\\Workspaces\\mgfsk_viterbi", 'mlse_vs_mmse_ber.png')
    shutil.copy(save_path, ws_path)
    
    out.append("\n" + "=" * 60)
    out.append("  仿真结果总结")
    out.append("=" * 60)
    out.append(f"  图片已保存: {save_path}")
    out.append("\n  关键观察：")
    out.append("  1. MLSE 在 ISI 信道中始终优于 MMSE 线性均衡")
    out.append("  2. 高 SNR 时，MLSE 的优势更明显（无噪声放大）")
    out.append("  3. MMSE 在强 ISI 下存在误差 floor")
    
    return {
        "output_text": "\n".join(out),
        "plot_path": save_path,
        "ws_plot_path": ws_path,
        "ber_mls": [float(x) for x in ber_mls],
        "ber_mmse": [float(x) for x in ber_mmse],
        "snr_range": [int(x) for x in SNR_DB_RANGE]
    }
