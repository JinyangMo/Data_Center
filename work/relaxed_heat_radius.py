import csv,json,math,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent;SRC=ROOT/'outputs'/'heat_all';OLD=ROOT/'outputs'/'heat_radius_only';OUT=ROOT/'outputs'/'heat_radius_relaxed';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader(open(SRC/'all_4590_records_heat_summary.csv',encoding='utf-8-sig')))
old={r['编号']:r for r in csv.DictReader(open(OLD/'每个数据中心_heat与粗略半径.csv',encoding='utf-8-sig'))}
profiles={}
for p in csv.DictReader(open(SRC/'period_temperature_profiles.csv',encoding='utf-8-sig')):profiles[(p['period_key'],int(p['inner_m']))]=p
RINGS=[(0,200),(200,400),(400,600),(600,800),(800,1000),(1000,1500),(1500,2000),(2000,2500),(2500,3000),(3000,4000),(4000,5000),(5000,6000)]
out=[];cross=collections.Counter();newcounts=collections.Counter();obs=collections.Counter()
for r in rows:
 prior=old[r['id']];obs[int(r.get('reference_scene_count') or 0)]+=1
 ps=[profiles.get((r.get('reference_period_key'),lo),{}) for lo,hi in RINGS]
 vals=[float(p['median_anomaly_c']) if p.get('median_anomaly_c') not in [None,''] and int(p.get('n_dates',0))>=1 else None for p in ps]
 core=[v for v in vals[:2] if v is not None];centre=sum(core)/len(core) if core else None
 band=None;reason='';method='';radius=None;lower=None;upper=None;rewarm=None;missing=None
 if prior['粗略heat半径_km_按0.5km取整']:
  radius=float(prior['粗略heat半径_km_按0.5km取整']);lower=float(prior['粗略距离带下限_km']);upper=float(prior['粗略距离带上限_km']);method='保留原规则结果';reason='已给出半径'
 elif not core:reason='近区没有有效温度，不能据此估算'
 elif centre<=1:reason='近区温差未超过1°C，未检出符合定义的热中心'
 else:
  for i in range(1,len(vals)-1):
   if vals[i] is not None and vals[i+1] is not None and abs(vals[i])<=1 and abs(vals[i+1])<=1:
    band=RINGS[i];lower=math.floor(band[0]/1000);upper=math.ceil(band[1]/1000);radius=(lower+upper)/2;method='一次观测也可用；首次连续两段回归；1km区间粗估';reason='已给出半径';rewarm=any(v is not None and v>1 for v in vals[i+2:]);missing=any(v is None for v in vals[i+2:]);break
  if band is None:reason='已有近区热中心，但现有曲线未找到连续两段回归，不能定半径'
 if r['coordinate_valid']=='False':reason='坐标缺失或无效'
 if radius is not None:
  assert lower<=radius<=upper and radius>0
  if band is not None:assert vals[i] is not None and vals[i+1] is not None
 newcounts[method or reason]+=1;cross[(prior['结果状态'],method or reason)]+=1
 out.append({'编号':r['id'],'data_center':r['site_name'],'纬度':r['latitude'],'经度':r['longitude'],'原规则半径_km':prior['粗略heat半径_km_按0.5km取整'],'建议使用半径_km':radius,'估计距离区间下限_km':lower,'估计距离区间上限_km':upper,'计算方式':method,'结果说明':reason,'使用宽松规则':band is not None,'宽松规则首个回归带下限_m':band[0] if band else None,'宽松规则首个回归带上限_m':band[1] if band else None,'后续是否再高于背景1C':rewarm,'后续是否仍有缺测':missing,'近区可用距离带数':len(core),'近区相对背景温差_C':round(centre,1) if centre is not None else None,'合格影像数':r.get('reference_scene_count') or 0,'位置是否近似':r['approximate_location']})
with open(OUT/'每个数据中心_宽松半径估算.csv','w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
stats={'all_records':len(rows),'original_numeric':522,'relaxed_added':sum(r['使用宽松规则'] for r in out),'total_numeric':sum(r['建议使用半径_km'] is not None for r in out),'relaxed_with_far_rewarming':sum(r['后续是否再高于背景1C'] is True for r in out),'relaxed_with_later_missing_data':sum(r['后续是否仍有缺测'] is True for r in out),'usable_scene_counts':dict(obs),'outcome_counts':dict(newcounts),'old_to_new':[{'old':a,'new':b,'n':n} for (a,b),n in cross.items()]}
assert len(out)==4590 and len({r['编号'] for r in out})==4590
assert stats['total_numeric']==stats['original_numeric']+stats['relaxed_added']
(OUT/'计算汇总.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'宽松估算口径.txt').write_text('只使用已有卫星地表温度，不使用建成年份，也不做DID。\n保留原规则能够识别的522条数值。其余记录：只要至少一个0—400米近区环带有一次合格观测且近区平均相对背景温差>1°C，就向外寻找首次连续两个环带都进入背景±1°C的位置。缺测不能当作回归。\n新估计按1公里距离区间表达，代表半径取该粗区间中点。距离区间是人为取整区间，不是置信区间或精度保证。该值是局部首次回归距离，不保证整个更远区域都恢复背景；后续再升温或缺测已另列字段。\n原规则值和宽松估算值的可靠程度不同；宽松结果可能受到云残留、邻近热源或地表类型变化影响。没有近区温度、未检出超过1°C的中心或没有观察到回归时，不生成数字，不用零或6公里替代。\n',encoding='utf-8')
print(json.dumps(stats,ensure_ascii=True))
