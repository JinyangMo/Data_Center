import csv
import json
import math
import statistics
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parent


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def describe(values):
    return {
        "数量": len(values), "平均半径_km": statistics.mean(values),
        "中位数_km": statistics.median(values),
        "最小值_km": min(values), "最大值_km": max(values),
    }


def analyze_distributions():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    original = read_csv(ROOT / "data/US Sites 2026-April (1).csv")
    rows = read_csv(ROOT / "outputs/heat_radius_complete/全部4590条_半径与来源.csv")
    strict_rows = read_csv(ROOT / "outputs/heat_radius_only/每个数据中心_heat与粗略半径.csv")
    original_ids = [r["id"] for r in original]
    assert len(set(original_ids)) == len(original_ids), "原始ID重复"
    assert [r["编号"] for r in rows] == original_ids, "全量半径与原始记录顺序或ID不一致"
    assert [r["编号"] for r in strict_rows] == original_ids, "严格结果与原始记录不一致"

    labels = ["卫星计算_原严格规则", "卫星计算_宽松首次回归", "模型估算_中位数插补"]
    groups = {label: [] for label in labels}
    full = []
    strict_by_id = {}
    for r in strict_rows:
        value = r["粗略heat半径_km_按0.5km取整"]
        if value != "":
            strict_by_id[r["编号"]] = float(value)
    for r in rows:
        value = float(r["最终半径_km"])
        assert math.isfinite(value) and value > 0
        assert r["半径来源"] in groups, "出现新的半径来源，请更新分类"
        groups[r["半径来源"]].append(value)
        full.append(value)
        if r["半径来源"] == labels[0]:
            assert strict_by_id[r["编号"]] == value, "严格半径发生变化"
    strict = list(strict_by_id.values())
    assert len(strict) == len(groups[labels[0]])
    assert sum(map(len, groups.values())) == len(full)
    # 当前结果均为0.5公里档。直接统计离散档位，避免直方图边界混淆。
    assert all(abs(v * 2 - round(v * 2)) < 1e-8 for v in full)
    centers = [i / 2 for i in range(1, math.ceil(max(full) * 2) + 1)]
    counters = {k: Counter(v) for k, v in groups.items()}
    frequencies = []
    for radius in centers:
        frequencies.append({
            "半径档位_km": radius,
            "全量数量": sum(c[radius] for c in counters.values()),
            "全量占比": full.count(radius) / len(full),
            "严格规则数量": counters[labels[0]][radius],
            "严格规则占比": strict.count(radius) / len(strict),
            "宽松新增数量": counters[labels[1]][radius],
            "模型插补数量": counters[labels[2]][radius],
        })
    assert sum(r["全量数量"] for r in frequencies) == len(original)
    assert sum(r["严格规则数量"] for r in frequencies) == len(strict)
    summary = {"全部最终半径": describe(full), "严格规则半径": describe(strict)}
    satellite = groups[labels[0]] + groups[labels[1]]
    summary["仅卫星半径"] = describe(satellite)
    summary["模型插补半径"] = describe(groups[labels[2]])
    summary["各来源数量"] = {k: len(v) for k, v in groups.items()}
    summary["说明"] = "按原始记录计数，可能共享园区或坐标。全量含中位数插补，不能当作全量实测热影响分布。严格平均使用表中粗略取整后的半径；未识别记录没有按零纳入。"

    out = ROOT / "outputs/radius_distribution"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "半径分布频数.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(frequencies[0]))
        writer.writeheader(); writer.writerows(frequencies)
    (out / "半径分布汇总.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    available = {f.name for f in font_manager.fontManager.ttflist}
    plt.rcParams["font.sans-serif"] = [f for f in ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"] if f in available]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    colors = ["#315A80", "#8AB5CE", "#E1B45B"]
    bottom = [0] * len(centers)
    for label, color, legend in zip(labels, colors, ["严格规则卫星结果", "宽松规则新增卫星结果", "模型插补"]):
        counts = [counters[label][x] for x in centers]
        axes[0].bar(centers, counts, bottom=bottom, width=.4, color=color, label=legend)
        bottom = [a + b for a, b in zip(bottom, counts)]
    axes[0].set_title(f"全部最终半径分布  n = {len(full):,}", fontsize=13, pad=13)
    strict_counts = [counters[labels[0]][x] for x in centers]
    axes[1].bar(centers, strict_counts, width=.4, color=colors[0], label="严格规则卫星结果")
    axes[1].set_title(f"严格规则半径分布  n = {len(strict):,}", fontsize=13, pad=13)
    for ax, values in zip(axes, [full, strict]):
        mean = statistics.mean(values)
        ax.axvline(mean, color="#AF3C39", ls="--", lw=1.7, label=f"平均半径 {mean:.2f} km")
        ax.set(xlabel="报告半径档位（km）", ylabel="记录数", xticks=centers)
        ax.set_ylim(0, max(bottom if ax is axes[0] else strict_counts) * 1.28)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.18);ax.set_axisbelow(True)
        ax.legend(loc="upper right", fontsize=9, frameon=False)
    for ax, counts in zip(axes, [bottom, strict_counts]):
        for x, n in zip(centers, counts):
            if n: ax.text(x, n, str(n), ha="center", va="bottom", fontsize=8)
    fig.text(.5, .025, "全量含模型插补：3,024条统一填入2.5 km；严格均值只使用通过严格规则的记录。两图纵轴尺度不同。", ha="center", fontsize=10)
    fig.subplots_adjust(left=.065, right=.98, bottom=.17, top=.87, wspace=.24)
    fig.savefig(out / "全部与严格规则半径分布.png", dpi=180)
    plt.close(fig)
    for name in ["全部最终半径", "严格规则半径", "仅卫星半径"]:
        result = summary[name]
        print(f"{name}：{result['数量']}条，平均 {result['平均半径_km']:.4f} km，中位数 {result['中位数_km']:.2f} km")
    print("输出目录：", out)
    return summary


if __name__ == "__main__":
    analyze_distributions()
