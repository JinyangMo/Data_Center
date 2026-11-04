import sys,os
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'heat_libs'))
import csv,json,math,collections,datetime
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'outputs'/'heat_all';CACHE=ROOT/'work'/'heat_all_cache'
source=list(csv.DictReader(open(ROOT/'outputs'/'data_centers_cleaned.csv',encoding='utf-8-sig')))
summary=list(csv.DictReader(open(OUT/'all_4590_records_heat_summary.csv',encoding='utf-8-sig')))
plan=json.loads((CACHE/'plan.json').read_text(encoding='utf-8'))
assert len(source)==len(summary)==len(plan['sites'])==4590
with open(OUT/'全部4590条_核心结果.csv',encoding='utf-8-sig') as f:core=list(csv.DictReader(f))
assert len(core)==4590 and [r['编号'] for r in core]==[r['id'] for r in source]
with open(OUT/'分析可用性分组.csv',encoding='utf-8-sig') as f:groups=list(csv.DictReader(f))
assert sum(int(r['原始记录数']) for r in groups)==4590
assert len({r['id'] for r in summary})==4590
for a,b,s in zip(source,summary,plan['sites']):
 assert a['id']==b['id']==s['id']
 assert all(a[k]==b[k] for k in a),'A source value changed'
 assert b['all_data_analysis_included']=='True'
 assert b['reference_status']!='pending'
 assert not b['formal_effect_radius_m']
 if s['coordinate_valid']:assert 'reference' in s['period_keys']
 else:assert not s['period_keys']
 if s['year_built'] is None:assert set(s['period_keys'])<= {'reference'}
 if s['year_built'] and s['year_built']>2026:assert set(s['period_keys'])<= {'reference'}
 for label,k in s['period_keys'].items():
  y=plan['jobs'][k]['year']
  if label=='reference':assert y==2026
  if label=='before':assert y==s['year_built']-2 and y>=1982
  if label=='after':assert y>s['year_built'] and y<=2026
stats=json.loads((OUT/'run_summary.json').read_text(encoding='utf-8'))
assert stats['completed_unique_coordinate_years']==len(plan['jobs'])==7293
periods={}
with open(OUT/'period_scene_results.jsonl',encoding='utf-8') as f:
 for line in f:
  p=json.loads(line);periods[p['key']]=p
assert set(periods)==set(plan['jobs'])
scene_count=0
for k,p in periods.items():
 job=plan['jobs'][k]
 assert p['year']==job['year']
 assert abs(p['latitude']-job['latitude'])<1e-8 and abs(p['longitude']-job['longitude'])<1e-8
 dates=[a['date'][:10] for a in p['accepted']];assert len(dates)==len(set(dates))<=2
 if p['status']=='ok':assert len(dates)==2
 for a in p['accepted']:
  d=datetime.date.fromisoformat(a['date'][:10]);assert d.year==p['year'] and 6<=d.month<=8
  assert a['scene_id'].split('_')[3]==d.strftime('%Y%m%d')
  assert len(a['bands'])==12
  assert a['inner_land_fraction']>=.25 and a['background_land_fraction']>=.1
  assert 0<=a['median_st_qa_k']<=5
  for b in a['bands']:
   assert 0<=b['coverage']<=1 and 0<=b['land_fraction']<=1
   if b['lst_c'] is not None:
    assert b['coverage']>=.7 and b['land_fraction']>=.25
    assert math.isfinite(b['lst_c']) and b['lst_c']>-273.15
    assert abs(b['lst_c']-b['background_c']-b['anomaly_c'])<1e-8
   else:assert b['anomaly_c'] is None
  scene_count+=1
with open(OUT/'ring_observations.csv',encoding='utf-8-sig') as f:observations=sum(1 for _ in csv.DictReader(f))
with open(OUT/'period_temperature_profiles.csv',encoding='utf-8-sig') as f:profiles=sum(1 for _ in csv.DictReader(f))
with open(OUT/'site_period_links.csv',encoding='utf-8-sig') as f:links=list(csv.DictReader(f))
assert observations==12*scene_count==stats['ring_observations']
assert profiles==12*len(plan['jobs'])
assert len(links)==sum(len(s['period_keys']) for s in plan['sites'])
assert {r['id'] for r in links}=={s['id'] for s in plan['sites'] if s['coordinate_valid']}
geo=json.loads((OUT/'all_valid_locations.geojson').read_text(encoding='utf-8'))
assert len(geo['features'])==4589 and {f['properties']['id'] for f in geo['features']}=={s['id'] for s in plan['sites'] if s['coordinate_valid']}
note=f'''数据核对完成
原始记录与全量结果均为4,590行，逐字段核对所有源值保持不变，id没有丢失或重复。
4,589条有效坐标全部建立参考观测任务；缺失年份、近似位置和所有状态均保留。
7,293个计划坐标—年份组合全部有明确处理结果，无pending。
核对了影像日期与场景标识、历史窗口、卫星档案范围、水陆覆盖、温度质量、温差算术、表间关联和GIS点数量。
共{scene_count:,}次不同坐标—年份的影像观测、{observations:,}条逐景环带观测和{profiles:,}条聚合环带记录。
同一卫星场景可能覆盖多个坐标，原始记录也可能共享坐标；这些数量均不代表独立热源或独立统计样本。
正式因果热影响半径全部保持缺失，没有把无法识别写成零。
本核对确认输出与计算规则一致，不证明投运年份、园区边界、匹配对照或因果解释正确。
'''
(OUT/'数据核对结果.txt').write_text(note,encoding='utf-8');print(note)
