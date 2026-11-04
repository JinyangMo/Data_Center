import csv,json,math,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent;SRC=ROOT/'outputs'/'heat_all';OUT=ROOT/'outputs'/'heat_radius_only';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader(open(SRC/'all_4590_records_heat_summary.csv',encoding='utf-8-sig')))
profiles={}
for p in csv.DictReader(open(SRC/'period_temperature_profiles.csv',encoding='utf-8-sig')):profiles[(p['period_key'],int(p['inner_m']))]=p
def numeric(r,k):return float(r[k]) if r.get(k) not in [None,''] else None
def rounded(v):return round(v,1) if v is not None else None
status_cn={'descriptive_recovery_found':'已识别粗略热区半径','insufficient_dates':'合格观测次数不足','insufficient_ring_dates':'部分距离带的有效观测不足','central_anomaly_not_consistently_above_threshold':'近区未形成持续高于背景1°C的热中心','not_stably_recovered_within_tested_range':'未满足稳定回归条件，半径未定'}
out=[];heat=[]
for r in rows:
 status=r.get('reference_recovery_1c_status','');lo=numeric(r,'reference_recovery_1c_recovery_band_inner_m');hi=numeric(r,'reference_recovery_1c_recovery_band_outer_m');rad=None;lower=None;upper=None
 if status=='descriptive_recovery_found':
  assert lo is not None and hi is not None
  rad=math.floor(((lo+hi)/2000)/.5+.5)*.5;lower=math.floor(lo/500)*.5;upper=math.ceil(hi/500)*.5
 reason=status_cn.get(status,'坐标缺失或无效' if r['coordinate_valid']=='False' else '温度观测不可用')
 out.append({'编号':r['id'],'data_center':r['site_name'],'纬度':r['latitude'],'经度':r['longitude'],'位置是否近似':r['approximate_location'],'近区地表温度_C':rounded(numeric(r,'reference_near_lst_c')),'背景地表温度_C':rounded(numeric(r,'reference_background_c')),'近区相对背景温差_C':rounded(numeric(r,'reference_near_anomaly_c')),'粗略heat半径_km_按0.5km取整':rad,'粗略距离带下限_km':lower,'粗略距离带上限_km':upper,'原始检出环带下限_m':lo,'原始检出环带上限_m':hi,'结果状态':reason,'影像采集日期':r.get('reference_dates','')})
 for a,b in [(0,200),(200,400),(400,600),(600,800),(800,1000),(1000,1500),(1500,2000),(2000,2500),(2500,3000),(3000,4000),(4000,5000),(5000,6000)]:
  p=profiles.get((r.get('reference_period_key'),a),{})
  heat.append({'编号':r['id'],'data_center':r['site_name'],'距离带内缘_m':a,'距离带外缘_m':b,'地表温度_C':rounded(numeric(p,'median_lst_c')),'相对背景温差_C':rounded(numeric(p,'median_anomaly_c')),'有效观测次数':int(p['n_dates']) if p else 0})
def write(name,data):
 with open(OUT/name,'w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
write('每个数据中心_heat与粗略半径.csv',out);write('每个数据中心_周边heat数据.csv',heat)
assert len(out)==4590 and len({r['编号'] for r in out})==4590 and len(heat)==4590*12
assert sum(r['粗略heat半径_km_按0.5km取整'] is not None for r in out)==522
assert all(r['粗略距离带下限_km']<=r['粗略heat半径_km_按0.5km取整']<=r['粗略距离带上限_km'] for r in out if r['粗略heat半径_km_按0.5km取整'] is not None)
note='只使用同一批2026年6—8月卫星地表温度数据，不清洗或使用建成年份，不进行DID或前后对比。\n粗略热区半径的操作定义：近区两个环带均高于6—8公里背景1°C，向外首次进入背景±1°C范围，且至少三个连续环带及所有更远已测环带都保持在范围内；各环带要求两次有效观测。半径取检出环带中点，按0.5公里取整，同时保留距离带。取整粒度不是测量精度保证。\n距离从数据中心提供的坐标点算起，搜索到6公里。未识别不能解释为零半径或半径大于6公里。空白温度表示缺测；有效观测次数为零是次数计数。\n这里的“heat辐射半径”仅表示卫星地表热区的粗略范围，未区分附近其他热源，不是对数据中心废热传播边界的确认。\n'
(OUT/'计算口径.txt').write_text(note,encoding='utf-8')
print(json.dumps({'rows':len(out),'numeric_radii':sum(r['粗略heat半径_km_按0.5km取整'] is not None for r in out),'statuses':dict(collections.Counter(r['结果状态'] for r in out))},ensure_ascii=True))
