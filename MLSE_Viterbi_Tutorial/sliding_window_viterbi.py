import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, shutil

# ============================================================
# 滑动窗口 Viterbi 实时检测演示
# ============================================================

def main(ctx):
    # 信道参数
    h = np.array([0.8, 0.5, 0.3], dtype=complex)
    L = len(h)
    M = 2
    sym_map = np.array([1, -1], dtype=complex)
    n_states = M ** (L - 1)
    
    # 构建 Trellis
    transitions = {}
    for s in range(n_states):
        bits = [(s >> (L - 2 - i)) & 1 for i in range(L - 1)]
        syms = sym_map[bits]
        for inp in [0, 1]:
            inp_sym = sym_map[inp]
            full = np.concatenate([[inp_sym], syms])
            expected = np.dot(h, full)
            next_s = (inp << (L - 2)) | (s >> 1)
            transitions[(s, inp)] = (next_s, expected)
    
    # 判决深度（典型值：5L）
    DECISION_DEPTH = 5 * L
    
    # 发送序列（模拟实时到达的数据）
    np.random.seed(42)
    N_TOTAL = 100
    tx_bits = np.random.randint(0, 2, N_TOTAL)
    tx_symbols = np.array([1 - 2*b for b in tx_bits], dtype=complex)
    
    # 通过信道加噪声
    r_clean = np.convolve(tx_symbols, h, mode='full')[:N_TOTAL]
    signal_power = np.mean(np.abs(r_clean) ** 2)
    snr_db = 8
    snr_linear = 10 ** (snr_db / 10)
    noise_power = signal_power / snr_linear
    noise = np.sqrt(noise_power / 2) * (np.random.randn(N_TOTAL) + 1j * np.random.randn(N_TOTAL))
    r = r_clean + noise
    
    # 前导：假设前 L-1 个符号是已知的 0（+1）
    prefix = np.zeros(L - 1, dtype=complex)
    r_extended = np.concatenate([prefix, r])
    N = len(r_extended)
    
    # 滑动窗口 Viterbi 实时检测
    INF = 1e18
    pm = np.full(n_states, INF)
    pm[0] = 0.0
    
    survivor_hist = []
    detected_bits = []
    output_times = []
    
    lines = []
    lines.append("=" * 70)
    lines.append("  Sliding-Window Viterbi (Real-Time MLSE)")
    lines.append("=" * 70)
    lines.append(f"\nChannel: h = {h}")
    lines.append(f"Channel memory L: {L}")
    lines.append(f"Decision depth D = 5L = {DECISION_DEPTH}")
    lines.append(f"Total symbols: {N_TOTAL}")
    lines.append(f"SNR: {snr_db} dB")
    lines.append(f"\nThe system has a fixed latency of D = {DECISION_DEPTH} symbol periods.")
    lines.append(f"Symbol at time t is output at time t + D.")
    
    lines.append(f"\n{'-' * 70}")
    lines.append(f"  Real-time processing log (first 20 steps)")
    lines.append(f"{'-' * 70}")
    
    for t in range(N):
        npm = np.full(n_states, INF)
        surv = np.zeros(n_states, dtype=int) - 1
        
        for s in range(n_states):
            if pm[s] >= INF:
                continue
            for inp in [0, 1]:
                ns, exp = transitions[(s, inp)]
                bm = np.abs(r_extended[t] - exp) ** 2
                cand = pm[s] + bm
                if cand < npm[ns]:
                    npm[ns] = cand
                    surv[ns] = s
        
        pm = npm
        survivor_hist.append(surv)
        
        if t >= DECISION_DEPTH:
            best_state = int(np.argmin(pm))
            state = best_state
            output_bit = None
            for back in range(DECISION_DEPTH):
                idx = t - back
                ps = survivor_hist[idx][state]
                if ps < 0:
                    break
                if back == DECISION_DEPTH - 1:
                    output_bit = (state >> (L - 2)) & 1
                state = ps
            
            if output_bit is not None:
                output_time = t - DECISION_DEPTH - (L - 1)
                if output_time >= 0 and output_time < N_TOTAL:
                    detected_bits.append(output_bit)
                    output_times.append(output_time)
                    
                    if len(detected_bits) <= 20:
                        true_bit = tx_bits[output_time]
                        err = "X" if output_bit != true_bit else "OK"
                        lines.append(f"  t={t:3d}: output bit[{output_time:2d}] = {output_bit} (true={true_bit}) {err}")
        
        if len(survivor_hist) > DECISION_DEPTH:
            survivor_hist.pop(0)
    
    # 处理末尾
    lines.append(f"\n{'-' * 70}")
    lines.append(f"  End-of-frame traceback (last {DECISION_DEPTH} symbols)")
    lines.append(f"{'-' * 70}")
    
    final_state = int(np.argmin(pm))
    state = final_state
    for back in range(min(DECISION_DEPTH, N_TOTAL - len(detected_bits))):
        idx = len(survivor_hist) - 1 - back
        if idx < 0:
            break
        ps = survivor_hist[idx][state]
        if ps < 0:
            break
        output_bit = (state >> (L - 2)) & 1
        output_time = N_TOTAL - 1 - back
        if output_time >= 0:
            detected_bits.append(output_bit)
            output_times.append(output_time)
            true_bit = tx_bits[output_time]
            err = "X" if output_bit != true_bit else "OK"
            lines.append(f"  final traceback: bit[{output_time:2d}] = {output_bit} (true={true_bit}) {err}")
        state = ps
    
    detected_bits = np.array(detected_bits)
    output_times = np.array(output_times)
    sort_idx = np.argsort(output_times)
    detected_bits = detected_bits[sort_idx]
    output_times = output_times[sort_idx]
    
    valid_len = min(len(detected_bits), N_TOTAL)
    errors = np.sum(detected_bits[:valid_len] != tx_bits[:valid_len])
    ber = errors / valid_len if valid_len > 0 else 0
    
    lines.append(f"\n{'=' * 70}")
    lines.append(f"  Summary")
    lines.append(f"{'=' * 70}")
    lines.append(f"  Total symbols: {N_TOTAL}")
    lines.append(f"  Valid output symbols: {valid_len}")
    lines.append(f"  Errors: {errors}")
    lines.append(f"  BER: {ber:.4f}")
    lines.append(f"  System latency: {DECISION_DEPTH} symbol periods")
    lines.append(f"\n  Key insight:")
    lines.append(f"  The receiver outputs bit[t] at time t+{DECISION_DEPTH},")
    lines.append(f"  allowing real-time streaming without waiting for end-of-frame.")
    
    report_text = "\n".join(lines)
    
    ws_dir = "C:\\Users\\Administrator\\Documents\\Kimi\\Workspaces\\mgfsk_viterbi"
    report_path = os.path.join(ws_dir, "sliding_window_viterbi.txt")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    # 绘制 Trellis 汇聚示意图
    fig, ax = plt.subplots(figsize=(12, 6))
    np.random.seed(123)
    n_paths = 8
    D = DECISION_DEPTH
    
    all_paths = []
    for p in range(n_paths):
        path = [p % n_states]
        for t in range(D + 5):
            if t < 3:
                path.append(np.random.randint(0, n_states))
            else:
                path.append(2)
        all_paths.append(path)
    
    for p in range(n_paths):
        path = all_paths[p]
        t_vals = list(range(len(path)))
        ax.plot(t_vals[:4], path[:4], 'o-', alpha=0.3, color='gray', markersize=4, linewidth=1)
        ax.plot(t_vals[3:], path[3:], 'o-', alpha=0.8, color='#1f77b4', markersize=4, linewidth=2)
    
    ax.axvline(x=3, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Convergence point (~3L)')
    ax.axvspan(0, 3, alpha=0.1, color='yellow', label='Transient (unreliable)')
    ax.axvspan(3, D+5, alpha=0.1, color='green', label='Steady-state (reliable)')
    
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('State', fontsize=12)
    ax.set_yticks(range(n_states))
    ax.set_title('Trellis Path Convergence: Why Truncated Traceback Works', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plot_path = os.path.join(ctx["runDir"], 'trellis_convergence.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    ws_plot_path = os.path.join(ws_dir, 'trellis_convergence.png')
    shutil.copy(plot_path, ws_plot_path)
    
    return {
        "report_path": report_path,
        "plot_path": plot_path,
        "ws_plot_path": ws_plot_path,
        "ber": float(ber),
        "latency": DECISION_DEPTH
    }
