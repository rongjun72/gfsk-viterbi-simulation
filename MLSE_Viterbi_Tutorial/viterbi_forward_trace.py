import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================================
# 信道参数
# ============================================================
h = np.array([0.8, 0.5, 0.3], dtype=complex)
L = len(h)
M = 2
sym_map = np.array([1, -1], dtype=complex)
n_states = M ** (L - 1)

# ============================================================
# 构建 Trellis 转移表
# ============================================================
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

# ============================================================
# 构造一个具体的接收信号例子
# 假设发送比特: 0, 1, 0, 1 (对应符号 +1, -1, +1, -1)
# 从状态0 ([+,+]) 开始，无噪声时接收信号应为: 1.6, 0.0, 0.6, -1.0
# 我们加一些噪声，让 r = [1.5, 0.2, 0.7, -0.9]
# ============================================================
r = np.array([1.5, 0.2, 0.7, -0.9], dtype=complex)
N = len(r)

# ============================================================
# Viterbi 前向递推（详细记录）
# ============================================================
INF = 1e18

pm_history = np.zeros((N + 1, n_states)) + INF
pm_history[0, 0] = 0.0

survivor = np.zeros((N, n_states), dtype=int) - 1
survivor_input = np.zeros((N, n_states), dtype=int) - 1

detail_log = []

for t in range(N):
    log_entry = {
        'time': t,
        'received': r[t].real,
        'transitions': []
    }
    
    for s in range(n_states):
        if pm_history[t, s] >= INF:
            continue
        for inp in [0, 1]:
            ns, expected = transitions[(s, inp)]
            bm = np.abs(r[t] - expected) ** 2
            cand = pm_history[t, s] + bm
            
            log_entry['transitions'].append({
                'from_state': s,
                'input': inp,
                'to_state': ns,
                'expected': expected.real,
                'branch_metric': bm.real,
                'prev_pm': pm_history[t, s],
                'candidate': cand.real
            })
            
            if cand < pm_history[t + 1, ns]:
                pm_history[t + 1, ns] = cand
                survivor[t, ns] = s
                survivor_input[t, ns] = inp
    
    detail_log.append(log_entry)

# ============================================================
# 回溯
# ============================================================
final_state = int(np.argmin(pm_history[N, :]))
state = final_state
detected_bits = np.zeros(N, dtype=int)
for t in range(N - 1, -1, -1):
    ps = survivor[t, state]
    if ps < 0:
        break
    inp = survivor_input[t, state]
    detected_bits[t] = inp
    state = ps

# ============================================================
# 打印详细报告
# ============================================================
report = []
report.append("=" * 70)
report.append("  Viterbi 前向递推过程 逐步演示")
report.append("=" * 70)
report.append(f"\n信道: h = {h}")
report.append(f"接收信号: r = {[round(x.real, 3) for x in r]}")
report.append(f"时刻数: N = {N}")
report.append(f"状态数: {n_states} (M^(L-1) = 2^2 = 4)")
report.append(f"\n初始条件: t=0 时, 只有状态0的累计度量 = 0")
report.append(f"          状态1,2,3 的累计度量 = inf")

for t in range(N):
    log = detail_log[t]
    report.append(f"\n{'-' * 70}")
    report.append(f"  时刻 t = {t + 1} (接收信号 r[{t}] = {log['received']:.3f})")
    report.append(f"{'-' * 70}")
    
    for ns in range(n_states):
        candidates = [x for x in log['transitions'] if x['to_state'] == ns]
        if not candidates:
            continue
        
        report.append(f"\n    到达状态 {ns} 的候选路径:")
        for c in candidates:
            flag = "  <- 胜出(幸存者)" if c['candidate'] == min([x['candidate'] for x in candidates]) else ""
            report.append(f"      从状态 {c['from_state']} 输入 {'+1' if c['input'] == 0 else '-1':>3} | "
                         f"期望输出 = {c['expected']:>6.3f} | "
                         f"分支度量 = |{log['received']:.3f} - {c['expected']:.3f}|^2 = {c['branch_metric']:.4f} | "
                         f"累计 = {c['prev_pm']:.4f} + {c['branch_metric']:.4f} = {c['candidate']:.4f}{flag}")
        
        winner = min(candidates, key=lambda x: x['candidate'])
        report.append(f"    -> 状态 {ns} 的幸存者: 来自状态 {winner['from_state']}, 累计度量 = {winner['candidate']:.4f}")
    
    report.append(f"\n    >>> 时刻 {t+1} 结束后的累计度量: " + 
                 ", ".join([f"PM[{s}]={pm_history[t+1, s]:.4f}" if pm_history[t+1, s] < INF else f"PM[{s}]=inf" for s in range(n_states)]))

report.append(f"\n{'=' * 70}")
report.append("  终止与回溯")
report.append(f"{'=' * 70}")
report.append(f"\n  最终时刻 t={N}:")
report.append(f"    各状态累计度量: " + 
             ", ".join([f"PM[{s}]={pm_history[N, s]:.4f}" for s in range(n_states)]))
report.append(f"    最小度量状态: {final_state} (累计度量 = {pm_history[N, final_state]:.4f})")

report.append(f"\n  回溯过程:")
state = final_state
path_str = [f"状态{state}(t={N})"]
for t in range(N - 1, -1, -1):
    ps = survivor[t, state]
    inp = survivor_input[t, state]
    if ps < 0:
        break
    path_str.insert(0, f"状态{ps}(t={t}) --[{inp}]--> ")
    state = ps

report.append(f"    最优路径: " + "".join(path_str))
report.append(f"    检测比特: {[int(x) for x in detected_bits]}")
report.append(f"    (0=+1, 1=-1)")

report.append(f"\n{'=' * 70}")
report.append("  结果对比")
report.append(f"{'=' * 70}")
report.append(f"  发送假设: 0, 1, 0, 1 (+1, -1, +1, -1)")
report.append(f"  检测比特: {[int(x) for x in detected_bits]}")
report.append(f"  正确!" if np.all(detected_bits == np.array([0, 1, 0, 1])) else "  有误")

report_path = os.path.join("C:\\Users\\Administrator\\Documents\\Kimi\\Workspaces\\mgfsk_viterbi", "viterbi_forward_trace.txt")
with open(report_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(report))

# ============================================================
# 绘制 Trellis 路径图
# ============================================================
fig, ax = plt.subplots(figsize=(14, 8))

for t in range(N):
    for s in range(n_states):
        for inp in [0, 1]:
            ns, expected = transitions[(s, inp)]
            ax.plot([t, t + 1], [s, ns], 'k-', alpha=0.15, linewidth=0.5)

state = final_state
path_states = [state]
path_inputs = []
for t in range(N - 1, -1, -1):
    ps = survivor[t, state]
    inp = survivor_input[t, state]
    if ps < 0:
        break
    path_states.insert(0, ps)
    path_inputs.insert(0, inp)
    state = ps

for t in range(N):
    ax.plot([t, t + 1], [path_states[t], path_states[t + 1]], 
            'b-', linewidth=3, alpha=0.7, zorder=3)

for t in range(N + 1):
    for s in range(n_states):
        color = '#1f77b4' if path_states[t] == s else 'white'
        edge_color = '#1f77b4' if path_states[t] == s else 'gray'
        ax.scatter(t, s, s=300, c=color, edgecolors=edge_color, linewidths=2, zorder=4)
        ax.text(t, s, f'{s}', ha='center', va='center', fontsize=10, fontweight='bold',
               color='white' if path_states[t] == s else 'black', zorder=5)

for t in range(N):
    x = t + 0.5
    y = (path_states[t] + path_states[t + 1]) / 2
    ax.text(x, y + 0.25, f"{'+1' if path_inputs[t] == 0 else '-1'}", 
           ha='center', va='bottom', fontsize=11, color='blue', fontweight='bold', zorder=5)
    ax.text(x, y - 0.25, f"PM={pm_history[t+1, path_states[t+1]]:.2f}", 
           ha='center', va='top', fontsize=9, color='darkgreen', zorder=5)

for t in range(N):
    ax.text(t + 0.5, 3.5, f"r[{t}]={r[t].real:.2f}", ha='center', va='center', 
           fontsize=10, color='red', fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='red', alpha=0.8))

ax.set_xlim(-0.3, N + 0.3)
ax.set_ylim(-0.5, 3.8)
ax.set_xlabel('Time Step', fontsize=12)
ax.set_ylabel('State', fontsize=12)
ax.set_yticks(range(n_states))
ax.set_yticklabels([f'{s}' for s in range(n_states)])
ax.set_xticks(range(N + 1))
ax.set_xticklabels([f't={t}' for t in range(N + 1)])
ax.set_title('Viterbi Trellis: Survivor Path (Blue) vs All Possible Transitions (Gray)', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')
ax.axhline(y=-0.1, color='black', linewidth=0.5)

plt.tight_layout()
plot_path = os.path.join("C:\\Users\\Administrator\\Documents\\Kimi\\Workspaces\\mgfsk_viterbi", "viterbi_trellis_path.png")
plt.savefig(plot_path, dpi=150, bbox_inches='tight')
plt.close()

print("\n".join(report))
print(f"\n报告已保存: {report_path}")
print(f"路径图已保存: {plot_path}")

