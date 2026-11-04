import csv,json,sys,math,hashlib,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'work/heat_libs'))
import numpy as np
OUT=ROOT/'outputs/heat_radius_complete';OUT.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(open(p,encoding='utf-8-sig')))
raw=read(ROOT/'outputs/heat_all/all_4590_records_heat_summary.csv')
prev=read(ROOT/'outputs/heat_radius_relaxed/每个数据中心_宽松半径估算.csv')
assert [r['id'] for r in raw]==[r['编号'] for r in prev]
def num(v):
 try:
  x=float(v);return x if math.isfinite(x) else np.nan
 except:return np.nan
fields=['latitude','longitude','current_mw','total_mw','total_building_sqft','total_gross_sqft']
X=np.array([[num(r[f]) for f in fields] for r in raw])
X[:,2:]=np.where(X[:,2:]>=0,np.log1p(np.maximum(X[:,2:],0)),np.nan)
y=np.array([num(r['建议使用半径_km']) for r in prev]);known=np.flatnonzero(np.isfinite(y))
# Connected components prevent a campus OR shared coordinate crossing folds.
parent=list(range(len(raw)))
def find(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
seen={}
for i,r in enumerate(raw):
 keys=[]
 if r['campus_id']:keys.append('campus:'+r['campus_id'])
 if np.isfinite(X[i,:2]).all():keys.append('coord:'+','.join(f'{v:.4f}' for v in X[i,:2]))
 for key in keys:
  if key in seen:parent[find(i)]=find(seen[key])
  else:seen[key]=i
groups=np.array([find(i) for i in range(len(raw))])
fold=np.array([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%5 for g in groups])
def neighbors(train,test):
 z=X[train];med=np.nanmedian(z,axis=0);scale=np.nanpercentile(z,75,axis=0)-np.nanpercentile(z,25,axis=0);scale=np.where(scale>0,scale,1)
 a=np.clip((np.where(np.isfinite(X),X,med)-med)/scale,-5,5)
 # Geographic coordinates together get half the weight; four size features the other half.
 weights=np.array([.25,.25,.125,.125,.125,.125])
 d=((a[test,None,:]-a[None,train,:])**2*weights).sum(axis=2)
 d+=((~np.isfinite(X[test,None,:]))!=(~np.isfinite(X[None,train,:]))).dot(weights)*.5
 return np.argsort(d,axis=1)
ks=[15,30,60];preds={str(k):np.full(len(raw),np.nan) for k in ks};preds['median']=np.full(len(raw),np.nan);preds['mean']=np.full(len(raw),np.nan)
for f in range(5):
 tr=known[fold[known]!=f];te=known[fold[known]==f];order=neighbors(tr,te)
 for k in ks:preds[str(k)][te]=np.median(y[tr[order[:,:k]]],axis=1)
 preds['median'][te]=np.median(y[tr]);preds['mean'][te]=np.mean(y[tr])
mae={k:float(np.mean(abs(v[known]-y[known]))) for k,v in preds.items()}
best=min(mae,key=mae.get)
missing=np.flatnonzero(~np.isfinite(y));order=neighbors(known,missing)
result=y.copy();donors={}
for j,i in enumerate(missing):
 ids=known[order[j,:int(best)]] if best.isdigit() else known
 result[i]=np.median(y[ids]) if best!='mean' else np.mean(y[ids]);donors[i]=ids
rows=[]
for i,(r,p) in enumerate(zip(raw,prev)):
 sat=bool(p['建议使用半径_km']);ids=donors.get(i,[])
 model=round(float(result[i])*2)/2 if not sat else None
 source=('卫星计算_原严格规则' if p['原规则半径_km'] else '卫星计算_宽松首次回归') if sat else ('模型估算_中位数插补' if best=='median' else '模型估算_条件性插补')
 note=p['结果说明']
 if not sat:note+='；假设存在可检测热中心时的参照半径，不证明该站实际存在热影响'
 rows.append({'编号':r['id'],'data_center':r['site_name'],'最终半径_km':float(y[i]) if sat else model,'半径来源':source,'卫星计算半径_km':float(y[i]) if sat else None,'模型估算半径_km':model,'原卫星结果说明':p['结果说明'],'解释':note,'纬度':r['latitude'],'经度':r['longitude'],'位置是否近似':r['approximate_location'],'近区相对背景温差_C':p['近区相对背景温差_C'],'合格影像数':int(p['合格影像数']),'卫星后续再次偏热':p['后续是否再高于背景1C'],'模型参照记录数':len(ids) if not sat else None,'模型参照半径P10_km':float(np.quantile(y[ids],.1)) if not sat else None,'模型参照半径P90_km':float(np.quantile(y[ids],.9)) if not sat else None,'模型参照编号':','.join(raw[j]['id'] for j in ids) if not sat and best.isdigit() else ('全部1566条卫星半径记录' if not sat else '')})
assert len(rows)==4590 and len({r['编号'] for r in rows})==4590
assert all(math.isfinite(r['最终半径_km']) and r['最终半径_km']>0 for r in rows)
assert sum(r['卫星计算半径_km'] is not None for r in rows)==1566
assert all(r['最终半径_km']==y[i] for i,r in enumerate(rows) if np.isfinite(y[i]))
with open(OUT/'全部4590条_半径与来源.csv','w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
cv=[]
for i in known:cv.append({'编号':raw[i]['id'],'验证分组':int(groups[i]),'折':int(fold[i]),'卫星半径_km':float(y[i]),'留出预测_km':float(preds[best][i]),'绝对误差_km':float(abs(preds[best][i]-y[i]))})
with open(OUT/'模型交叉验证.csv','w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=list(cv[0]));w.writeheader();w.writerows(cv)
stats={'记录总数':len(rows),'卫星计算':1566,'其中原严格规则':522,'其中宽松规则':1044,'模型估算':len(missing),'全量混合平均_km':float(np.mean([r['最终半径_km'] for r in rows])),'卫星平均_km':float(np.mean(y[known])),'模型平均_km':float(np.mean([r['模型估算半径_km'] for r in rows if r['模型估算半径_km'] is not None])),'分组交叉验证MAE_km':mae,'选定模型':best,'验证误差中位数_km':float(np.median(abs(preds[best][known]-y[known]))),'验证误差P90_km':float(np.quantile(abs(preds[best][known]-y[known]),.9)),'训练记录数':len(known),'训练独立组数':len(set(groups[known].tolist())),'模型填补原始原因':dict(collections.Counter(prev[i]['结果说明'] for i in missing))}
(OUT/'模型与结果汇总.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2),encoding='utf-8')
notes=[['项目','说明'],['卫星时期','2026年6—8月。仅使用现有结果，不处理年份、不做DID。'],['覆盖','4590条记录全部有数值。1566条保留卫星计算，3024条模型插补。'],['卫星半径','522条原严格规则，1044条宽松首次局部回归。并非独立实测废热扩散边界。'],['预测目标','在可检测热中心且可识别回归距离的参照样本中，预测粗略首次回归半径。'],['预测变量','仅纬度、经度、当前功率、总功率、建筑面积、总建筑面积。不用年份或建成前后结果。'],['方法','比较15、30、60近邻中位数与全局中位数、均值。按分组五折留出MAE最小选取。选定：'+best],['距离规则','规模取log(1+x)，训练集IQR标准化并截断±5。坐标合计权重0.5，四项规模合计0.5。缺失以训练中位数填充，缺失状态不一致加0.5倍相应权重。'],['防止验证泄漏','同一campus_id或四位小数坐标相同的记录连通分组，整组划入同一折。'],['重复站点','源文件按记录保留。相同坐标/园区可能重复，训练和全量平均仍按记录加权。'],['模型报告精度','模型半径取最近0.5公里，仅为粗略表达，不表示精度达到0.5公里。'],['参照P10/P90','参照记录半径分布的10%和90%分位数，不是预测置信区间。'],['未检出热中心','仍保留原卫星标记。模型正半径是假设性参照，不证明该站确实存在热影响，也不代表原观测为零。'],['无回归边界','未把900条简单设为大于6公里。非单调、缺测、低于背景或邻近热源均可能导致原规则不通过。模型并未解决这些成因。'],['缺坐标记录','只使用其可用规模信息或训练中位数，仍有数值但缺乏位置依据。'],['适用范围限制','训练仅来自能识别半径的1566条，属于有选择的样本。模型对未检出热中心/无回归/缺测记录的误差无法由现有数据验证。'],['验证解释','验证仅测量模型对已有卫星半径的复现误差，不代表真实物理边界误差。模型复杂度选择使用同一五折，结果属于内部验证。'],['全量平均解释','混合平均含3024条模型插补，不能称为全量卫星实测平均或数据中心造成的平均热影响。'],['原始来源','US Sites 2026-April (1).csv；已有Landsat Collection 2 Level-2地表温度提取结果。'],['卫星产品说明','https://www.usgs.gov/landsat-missions/landsat-collection-2-surface-temperature']]
(OUT/'方法说明.txt').write_text('\n\n'.join(a+'：'+b for a,b in notes[1:]),encoding='utf-8')
(ROOT/'work/radii_workbook_data.json').write_text(json.dumps({'rows':rows,'stats':stats,'notes':notes},ensure_ascii=False),encoding='utf-8')
print(json.dumps(stats,ensure_ascii=False))
