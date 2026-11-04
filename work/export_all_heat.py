import sys,os
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'heat_libs'))
import csv,json,math,collections,datetime
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'outputs'/'heat_all';CACHE=ROOT/'work'/'heat_all_cache'
plan=json.loads((CACHE/'plan.json').read_text(encoding='utf-8'));rows=list(csv.DictReader(open(ROOT/'outputs'/'data_centers_cleaned.csv',encoding='utf-8-sig')))
look={r['id']:r for r in rows};periods={};profiles=[];observations=[];recoveries=[];pstats={}
def write(name,data):
 if not data:return
 keys=list(dict.fromkeys(k for row in data for k in row))
 with open(OUT/name,'w',newline='',encoding='utf-8-sig') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(data)
for k,job in plan['jobs'].items():
 path=CACHE/'periods'/(k+'.json')
 if not path.exists():continue
 try:p=json.loads(path.read_text(encoding='utf-8'))
 except json.JSONDecodeError:continue
 periods[k]=p;agg=[]
 for a in p['accepted']:
  for b in a['bands']:observations.append(dict(period_key=k,year=p['year'],longitude=p['longitude'],latitude=p['latitude'],scene_id=a['scene_id'],date=a['date'],platform=a.get('platform',a['scene_id'].split('_')[0]),median_inner_st_qa_k=a['median_st_qa_k'],**b))
 for i,(lo,hi) in enumerate([(0,200),(200,400),(400,600),(600,800),(800,1000),(1000,1500),(1500,2000),(2000,2500),(2500,3000),(3000,4000),(4000,5000),(5000,6000)]):
  bands=[a['bands'][i] for a in p['accepted'] if a['bands'][i]['lst_c'] is not None]
  z={'period_key':k,'year':p['year'],'inner_m':lo,'outer_m':hi,'n_dates':len(bands),'median_lst_c':float(np.median([b['lst_c'] for b in bands])) if bands else None,'median_anomaly_c':float(np.median([b['anomaly_c'] for b in bands])) if bands else None,'median_land_fraction':float(np.median([a['bands'][i]['land_fraction'] for a in p['accepted']])) if p['accepted'] else None}
  profiles.append(z);agg.append(z)
 near=agg[:2];pst={'status':p['status'],'year':p['year'],'scene_count':len(p['accepted']),'scene_ids':';'.join(a['scene_id'] for a in p['accepted']),'dates':';'.join(a['date'][:10] for a in p['accepted']),'platforms':';'.join(sorted({a.get('platform',a['scene_id'].split('_')[0]) for a in p['accepted']})),'median_st_qa_k':float(np.median([a['median_st_qa_k'] for a in p['accepted']])) if p['accepted'] else None,'median_day_of_year':float(np.median([datetime.date.fromisoformat(a['date'][:10]).timetuple().tm_yday for a in p['accepted']])) if p['accepted'] else None,'background_c':float(np.median([a['background_c'] for a in p['accepted']])) if p['accepted'] else None,'near_anomaly_c':float(np.mean([b['median_anomaly_c'] for b in near])) if all(b['median_anomaly_c'] is not None for b in near) else None,'near_lst_c':float(np.mean([b['median_lst_c'] for b in near])) if all(b['median_lst_c'] is not None for b in near) else None,'near_min_date_count':min(b['n_dates'] for b in near)}
 vals=[b['median_anomaly_c'] for b in agg[:9]];valid=[i for i,v in enumerate(vals) if v is not None]
 rho=None;decl=None
 if len(valid)>=6:
  x=np.arange(len(valid));v=np.array([vals[i] for i in valid]);_,inv,counts=np.unique(v,return_inverse=True,return_counts=True);ranks=(np.cumsum(counts)-counts+(counts-1)/2)[inv];rho=float(np.corrcoef(x,ranks)[0,1]) if np.std(ranks)>0 else 0.0
  pairs=[(vals[i],vals[i+1]) for i in range(8) if vals[i] is not None and vals[i+1] is not None];decl=sum(b<a for a,b in pairs)/len(pairs) if pairs else None
 pst['rank_correlation_0_3km']=rho;pst['adjacent_decline_fraction_0_3km']=decl;pst['decreasing_shape']=bool(p['status']=='ok' and rho is not None and rho<=-.7 and decl is not None and decl>=.625 and pst['near_anomaly_c'] is not None and pst['near_anomaly_c']>=1)
 for tol in [1,2,3]:
  rr={'period_key':k,'year':p['year'],'tolerance_c':tol,'status':'insufficient_dates','recovery_band_inner_m':None,'recovery_band_outer_m':None,'interpretation':'descriptive_point_based_recovery_not_causal_effect_radius'}
  if p['status']=='ok':
   if any(b['n_dates']<2 for b in agg):rr['status']='insufficient_ring_dates'
   elif min(b['median_anomaly_c'] for b in near)<=tol:rr['status']='central_anomaly_not_consistently_above_threshold'
   else:
    rr['status']='not_stably_recovered_within_tested_range'
    for i in range(1,10):
     if all(abs(b['median_anomaly_c'])<=tol for b in agg[i:]):rr.update(status='descriptive_recovery_found',recovery_band_inner_m=agg[i]['inner_m'],recovery_band_outer_m=agg[i]['outer_m']);break
  recoveries.append(rr)
  for key in ['status','recovery_band_inner_m','recovery_band_outer_m']:pst[f'recovery_{tol}c_{key}']=rr[key]
 pstats[k]=pst
summary=[];links=[];opids=[];oplat=[];oplon=[]
for s in plan['sites']:
 if s['coordinate_valid'] and 'Operational' in look[s['id']]['stage']:opids.append(s['id']);oplat.append(s['latitude']);oplon.append(s['longitude'])
oplat=np.radians(oplat);oplon=np.radians(oplon);opids=np.array(opids)
for s in plan['sites']:
 r={**look[s['id']],'all_data_analysis_included':True,'coordinate_valid':s['coordinate_valid'],'analysis_reference_year':plan['reference_year'],'formal_effect_radius_m':None,'formal_effect_radius_status':'not_estimated_no_verified_opening_or_matched_controls'}
 y=s['year_built'];r['year_group']='missing' if y is None else ('future_after_reference' if y>plan['reference_year'] else 'reference_year' if y==plan['reference_year'] else 'historical')
 r['current_stage_contains_operational']='Operational' in r['stage'];r['historical_pair_status']='not_planned'
 for label in ['reference','before','after']:
  k=s['period_keys'].get(label)
  if k:
   links.append({'id':s['id'],'site_name':s['site_name'],'period_role':label,'period_key':k,'year':plan['jobs'][k]['year']})
   r[label+'_period_key']=k
   if k in pstats:
    for name,value in pstats[k].items():r[label+'_'+name]=value
   else:r[label+'_status']='pending'
  else:r[label+'_status']=s['unavailable_reasons'].get(label,'not_planned')
 b=pstats.get(s['period_keys'].get('before'));a=pstats.get(s['period_keys'].get('after'))
 r['change_near_anomaly_c']=None;r['comparison_day_of_year_shift']=None;r['comparison_uses_different_platforms']=None
 if b and a:
  if b['status']=='ok' and a['status']=='ok':r['historical_pair_status']='two_dates_each_period' if b['near_min_date_count']>=2 and a['near_min_date_count']>=2 else 'incomplete_near_ring_coverage'
  else:r['historical_pair_status']='insufficient_dates'
  if b['near_anomaly_c'] is not None and a['near_anomaly_c'] is not None:
   r['change_near_anomaly_c']=a['near_anomaly_c']-b['near_anomaly_c']
  if b['median_day_of_year'] is not None and a['median_day_of_year'] is not None:r['comparison_day_of_year_shift']=a['median_day_of_year']-b['median_day_of_year']
  r['comparison_uses_different_platforms']=a['platforms']!=b['platforms']
 elif 'before' in s['period_keys'] and 'after' in s['period_keys']:r['historical_pair_status']='pending_periods'
 else:r['historical_pair_status']=';'.join(f'{k}:{v}' for k,v in s['unavailable_reasons'].items() if k in ['before','after']) or 'not_planned'
 r['interpretation']='descriptive_only_year_built_proxy_unverified_stage_not_past_operation_evidence'
 r['nearest_other_listed_operational_record_m']=None;r['other_listed_operational_records_within_8km']=None
 if s['coordinate_valid'] and len(opids):
  lat=math.radians(s['latitude']);lon=math.radians(s['longitude']);v=np.sin((oplat-lat)/2)**2+math.cos(lat)*np.cos(oplat)*np.sin((oplon-lon)/2)**2;dist=6371000*2*np.arcsin(np.sqrt(np.clip(v,0,1)));d=dist[opids!=s['id']]
  r['nearest_other_listed_operational_record_m']=float(np.min(d)) if len(d) else None;r['other_listed_operational_records_within_8km']=int(np.sum(d<8000))
 summary.append(r)
write('all_4590_records_heat_summary.csv',summary);write('site_period_links.csv',links);write('period_temperature_profiles.csv',profiles);write('ring_observations.csv',observations);write('apparent_recovery_distances.csv',recoveries)
scene_sources={}
for p in periods.values():
 for a in p['accepted']:
  entry=scene_sources.setdefault(a['scene_id'],{'scene_id':a['scene_id'],'acquisition_utc':a['date'],'platform':a.get('platform'),'stac_item_url':'https://planetarycomputer.microsoft.com/api/stac/v1/collections/landsat-c2-l2/items/'+a['scene_id'],'coordinate_year_uses':0,'extraction_date':'2026-11-01'})
  entry['coordinate_year_uses']+=1
write('source_scenes.csv',list(scene_sources.values()))
core=[]
status_cn={'ok':'取得两次合格观测','insufficient_scenes':'本轮不足两次合格观测','access_error':'数据访问失败','invalid_or_missing_coordinates':'坐标缺失或无效','pending':'处理中'}
for r in summary:
 pair=r['historical_pair_status']
 if not r['coordinate_valid']:pair_cn='坐标缺失或无效'
 elif pair=='two_dates_each_period':pair_cn='前后近区均有两次合格观测'
 elif pair=='incomplete_near_ring_coverage':pair_cn='影像可用，但近区环带覆盖不足'
 elif pair=='insufficient_dates':pair_cn='前后观测不足'
 elif pair=='pending_periods':pair_cn='历史时段处理中'
 elif r['year_group']=='missing':pair_cn='年份缺失，无法划分前后'
 elif r['year_group']=='future_after_reference':pair_cn='未来项目，尚无投运后观测'
 elif r['year_group']=='reference_year':pair_cn='2026年记录，尚无完整建成后年份'
 else:pair_cn='部分历史窗口早于卫星档案'
 core.append({'编号':r['id'],'名称':r['site_name'],'公司':r['company_name'],'源表项目状态':r['stage'],'建成年份_尚未核实为投运年份':r['year_built'],'纬度':r['latitude'],'经度':r['longitude'],'是否近似位置':r['approximate_location'],'2026温度数据状态':status_cn.get(r['reference_status'],r['reference_status']),'2026合格观测日期':r.get('reference_dates'),'2026近区地表温度_C':r.get('reference_near_lst_c'),'2026近区相对背景温差_C':r.get('reference_near_anomaly_c'),'前期年份':r.get('before_year'),'后期年份':r.get('after_year'),'历史对比状态':pair_cn,'前后相对温差变化_C_非因果':r.get('change_near_anomaly_c'),'2026表观回归距离带下限_m_1C阈值_非因果':r.get('reference_recovery_1c_recovery_band_inner_m'),'2026表观回归距离带上限_m_1C阈值_非因果':r.get('reference_recovery_1c_recovery_band_outer_m'),'表观回归距离状态':r.get('reference_recovery_1c_status'),'8公里内其他运营记录数量':r.get('other_listed_operational_records_within_8km')})
for c,r in zip(core,summary):
 c['2026近区最少有效日期数']=r.get('reference_near_min_date_count')
 for field in ['2026近区地表温度_C','2026近区相对背景温差_C','前后相对温差变化_C_非因果']:
  if c[field] is not None:c[field]=round(c[field],2)
 recovery_cn={'descriptive_recovery_found':'找到表观距离带（非因果半径）','insufficient_dates':'观测日期不足','insufficient_ring_dates':'部分环带观测不足','central_anomaly_not_consistently_above_threshold':'近区温差未持续高于阈值','not_stably_recovered_within_tested_range':'搜索范围内未稳定回归，不能推断影响超过6公里'}
 c['表观回归距离状态']=recovery_cn.get(c['表观回归距离状态'],c['表观回归距离状态'])
write('全部4590条_核心结果.csv',core)
availability={}
for r in summary:
 if not r['coordinate_valid']:label='坐标缺失或无效'
 elif r['year_group']=='missing':label='年份缺失：保留参考温度，无法划分前后'
 elif r['year_group']=='future_after_reference':label='2027年及以后：保留所在地参考温度，尚无未来投运后观测'
 elif r['year_group']=='reference_year':label='2026年：保留参考与可用前期温度，尚无完整建成后年份'
 elif r.get('before_period_key') and r.get('after_period_key'):label='具备建成年份前后两个卫星检索窗口'
 else:label='部分历史窗口早于卫星温度档案：保留参考温度与可用单侧观测'
 g=availability.setdefault(label,{'类别':label,'原始记录数':0,'参考年份取得两景影像数':0,'参考近区具有两次有效观测数':0,'前后近区均有两次有效观测数':0})
 g['原始记录数']+=1;g['参考年份取得两景影像数']+=r.get('reference_status')=='ok';g['参考近区具有两次有效观测数']+=r.get('reference_status')=='ok' and r.get('reference_near_min_date_count',0)>=2;g['前后近区均有两次有效观测数']+=r['historical_pair_status']=='two_dates_each_period'
write('分析可用性分组.csv',list(availability.values()))
with open(OUT/'period_scene_results.jsonl','w',encoding='utf-8') as f:
 for p in periods.values():f.write(json.dumps(p,ensure_ascii=False)+'\n')
features=[]
for s,r in zip(plan['sites'],summary):
 if not s['coordinate_valid']:continue
 props={k:r.get(k) for k in ['id','site_name','company_name','company_type','stage','year_built','approximate_location','reference_period_key','before_period_key','after_period_key','reference_status','reference_near_anomaly_c','reference_near_lst_c','reference_scene_count','reference_decreasing_shape','reference_recovery_1c_status','reference_recovery_1c_recovery_band_inner_m','reference_recovery_1c_recovery_band_outer_m','historical_pair_status','change_near_anomaly_c','formal_effect_radius_status']}
 features.append({'type':'Feature','geometry':{'type':'Point','coordinates':[s['longitude'],s['latitude']]},'properties':props})
(OUT/'all_valid_locations.geojson').write_text(json.dumps({'type':'FeatureCollection','features':features},ensure_ascii=False),encoding='utf-8')
stats={'source_records':len(rows),'summary_records':len(summary),'valid_coordinate_records':sum(s['coordinate_valid'] for s in plan['sites']),'reference_year':plan['reference_year'],'planned_unique_coordinate_years':len(plan['jobs']),'completed_unique_coordinate_years':len(periods),'reference_status_counts':dict(collections.Counter(r['reference_status'] for r in summary)),'reference_two_dates_complete_near_count':sum(r.get('reference_status')=='ok' and r.get('reference_near_min_date_count',0)>=2 for r in summary),'historical_pair_two_dates_each':sum(r['historical_pair_status']=='two_dates_each_period' for r in summary),'historical_pair_status_counts':dict(collections.Counter(r['historical_pair_status'] for r in summary)),'missing_year_reference_ok':sum(not r['year_built'] and r['reference_status']=='ok' for r in summary),'approximate_location_reference_ok':sum(r['approximate_location']=='true' and r['reference_status']=='ok' for r in summary),'reference_decreasing_shapes':sum(bool(r.get('reference_decreasing_shape')) for r in summary),'reference_recovery_at_1c':sum(r.get('reference_recovery_1c_status')=='descriptive_recovery_found' for r in summary),'unique_coordinate_year_scene_observations':sum(len(p['accepted']) for p in periods.values()),'unique_source_scene_ids':len({a['scene_id'] for p in periods.values() for a in p['accepted']}),'ring_observations':len(observations),'source_year_counts':dict(collections.Counter(r['year_group'] for r in summary)),'formal_effect_radii_estimated':0}
(OUT/'run_summary.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
counts=collections.Counter(r['reference_status'] for r in summary);labels={'ok':'Two usable scenes','insufficient_scenes':'Insufficient scenes','invalid_or_missing_coordinates':'Invalid / missing coordinates','access_error':'Access failed','pending':'In progress'};axes[0].barh([labels.get(k,k) for k in counts],list(counts.values()),color='#427ca5');axes[0].set(xlabel='Source records',title='All 4,590 records: reference-data status');axes[0].set_xlim(0,max(counts.values())*1.13)
for i,v in enumerate(counts.values()):axes[0].text(v+max(counts.values())*.015,i,str(v),va='center',fontsize=10)
vals=[r['reference_near_anomaly_c'] for r in summary if r.get('reference_status')=='ok' and r.get('reference_near_anomaly_c') is not None]
if vals:axes[1].hist(vals,bins=40,color='#d2704b');axes[1].axvline(0,color='#444',lw=1)
axes[1].set(xlabel='0–400 m temperature minus 6–8 km background (°C)',ylabel='Source records',title='2026 snapshots at all project stages')
fig.suptitle('2026 daytime surface-temperature screening\nUp to two scenes per location; descriptive temperature contrasts')
fig.savefig(OUT/'all_records_overview.png',dpi=150,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5.5),gridspec_kw={'width_ratios':[4,1]},layout='constrained')
for ax,west,title in [(axes[0],False,'Contiguous United States'),(axes[1],True,'Hawaii / western locations')]:
 group=[r for r in summary if r['coordinate_valid'] and (float(r['longitude'])<-130)==west]
 ax.scatter([float(r['longitude']) for r in group],[float(r['latitude']) for r in group],s=7,c='#c6c9ce',linewidths=0,label='No complete two-date near-area result')
 good=[r for r in group if r.get('reference_status')=='ok' and r.get('reference_near_anomaly_c') is not None]
 if good:sc=ax.scatter([float(r['longitude']) for r in good],[float(r['latitude']) for r in good],c=[r['reference_near_anomaly_c'] for r in good],s=9,vmin=-8,vmax=8,cmap='coolwarm',linewidths=0)
 ax.set(title=title,xlabel='Longitude',ylabel='Latitude');ax.grid(alpha=.12)
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
fig.colorbar(ScalarMappable(norm=Normalize(-8,8),cmap='coolwarm'),ax=axes.tolist(),label='0–400 m temperature minus 6–8 km background (°C)',extend='both',shrink=.8)
fig.suptitle('Source locations: 2026 relative surface temperature snapshots\nAll project stages; gray = incomplete near-area data; unmatched background')
fig.savefig(OUT/'all_locations_reference_temperature.png',dpi=170,bbox_inches='tight');plt.close(fig)
print(json.dumps(stats,ensure_ascii=False),flush=True)
