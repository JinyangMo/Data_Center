import argparse, csv, json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description="Data Center Heat Analysis")
    parser.add_argument("--recalculate", action="store_true", help="用已有卫星汇总重算半径，不下载影像")
    args = parser.parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    if args.recalculate:
        for script in ("heat_radius_only.py", "relaxed_heat_radius.py", "impute_all_radii.py"):
            subprocess.run([sys.executable, str(ROOT / "work" / script)], cwd=ROOT, check=True)
    original = ROOT / "data" / "US Sites 2026-April (1).csv"
    with original.open(encoding="utf-8-sig", newline="") as f:
        count = sum(1 for _ in csv.DictReader(f))
    stats_path = ROOT / "outputs/heat_radius_complete/模型与结果汇总.json"
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    print("Data Center Heat Analysis")
    print(f"原始数据：{count} 条")
    print(f"卫星计算：{stats['卫星计算']} 条；模型插补：{stats['模型估算']} 条")
    print(f"卫星平均：{stats['卫星平均_km']:.2f} km；混合平均：{stats['全量混合平均_km']:.2f} km")
    strict_file = ROOT / "outputs/heat_radius_only/每个数据中心_heat与粗略半径.csv"
    with strict_file.open(encoding="utf-8-sig", newline="") as stream:
        strict_radii = [float(row["粗略heat半径_km_按0.5km取整"])
                        for row in csv.DictReader(stream)
                        if row["粗略heat半径_km_按0.5km取整"] != ""]
    if strict_radii:
        print(f"严格规则平均半径：{sum(strict_radii) / len(strict_radii):.2f} km（{len(strict_radii)}条）")
    print("结果目录：", stats_path.parent)

if __name__ == "__main__":
    main()
