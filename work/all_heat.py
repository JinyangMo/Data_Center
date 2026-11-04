import sys,os
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'heat_libs'))
import csv,json,time,math,hashlib,threading,datetime,urllib.request,urllib.parse,urllib.error,concurrent.futures,collections
from pathlib import Path
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from pyproj import Transformer
ROOT=Path(__file__).resolve().parent.parent;OUT=ROOT/'outputs'/'heat_all';OUT.mkdir(exist_ok=True)
CACHE=ROOT/'work'/'heat_all_cache';CACHE.mkdir(exist_ok=True)
PER=CACHE/'periods';PER.mkdir(exist_ok=True)
BASE='https://planetarycomputer.microsoft.com/api/stac/v1';REF=2026
RINGS=[(0,200),(200,400),(400,600),(600,800),(800,1000),(1000,1500),(1500,2000),(2000,2500),(2500,3000),(3000,4000),(4000,5000),(5000,6000)]
AUTH_LOCK=threading.Lock();NET_LOCK=threading.Lock();NET_SEM=threading.Semaphore(8)
AUTH={'token':'','expiry':0,'checked':0};NEXT_NET=0
def get_json(url):
 global NEXT_NET
 for attempt in range(4):
  try:
   with NET_LOCK:
    delay=max(0,NEXT_NET-time.monotonic());NEXT_NET=max(time.monotonic(),NEXT_NET)+.12
   if delay:time.sleep(delay)
   with NET_SEM:
    req=urllib.request.Request(url,headers={'User-Agent':'DataCenterTemperatureResearch/1.0'})
    with urllib.request.urlopen(req,timeout=35) as resp:return json.load(resp)
  except urllib.error.HTTPError as e:
   if e.code not in [429,500,502,503,504] or attempt==3:raise
   try:delay=float(e.headers.get('Retry-After','5'))
   except ValueError:delay=5
   time.sleep(min(30,max(delay,2**attempt)))
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
def token():
 now=time.time()
 if AUTH['expiry']-now>300:return AUTH['token']
 with AUTH_LOCK:
  now=time.time()
  if AUTH['expiry']-now>300 or (AUTH['expiry']>now+30 and now-AUTH['checked']<120):return AUTH['token']
  j=get_json('https://planetarycomputer.microsoft.com/api/sas/v1/token/landsateuwest/landsat-c2')
  AUTH.update(token=j['token'],expiry=datetime.datetime.fromisoformat(j['msft:expiry'].replace('Z','+00:00')).timestamp(),checked=now)
 return AUTH['token']
def atomic(path,obj):
 tmp=path.with_suffix(f'.{os.getpid()}.{threading.get_ident()}.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False),encoding='utf-8');tmp.replace(path)
def key(lon,lat,year):return hashlib.sha256(f'{lon:.8f},{lat:.8f},{year},land_coverage_v2'.encode()).hexdigest()[:24]
def coords(r):
 try:
  lat=float(r['latitude']);lon=float(r['longitude'])
  return (lon,lat) if math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180 and (lat or lon) else None
 except (ValueError,TypeError):return None
def plan():
 rows=list(csv.DictReader(open(ROOT/'outputs'/'data_centers_cleaned.csv',encoding='utf-8-sig')));jobs={};sites=[]
 for r in rows:
  c=coords(r);y=int(r['year_built_normalized']) if r['year_built_normalized'] else None
  s={'id':r['id'],'site_name':r['site_name'],'year_built':y,'coordinate_valid':bool(c),'period_keys':{},'unavailable_reasons':{}}
  if c:s.update(longitude=c[0],latitude=c[1])
  years={'reference':REF}
  if y:
   if y<=REF:
    if y-2>=1982:years['before']=y-2
    else:s['unavailable_reasons']['before']='before_landsat_temperature_archive'
    if y<REF:
     post=y+2 if y+2<=REF else y+1
     if 1982<=post<=REF:years['after']=post
     else:s['unavailable_reasons']['after']='outside_available_archive'
    else:s['unavailable_reasons']['after']='no_complete_post_year_yet'
   else:
    s['unavailable_reasons']['before']='future_project_reference_is_baseline_only';s['unavailable_reasons']['after']='future_year_not_observed'
  else:
   s['unavailable_reasons']['before']='missing_year';s['unavailable_reasons']['after']='missing_year'
  if c:
   for label,year in years.items():
    k=key(*c,year);s['period_keys'][label]=k;jobs[k]={'key':k,'longitude':c[0],'latitude':c[1],'year':year}
  else:s['unavailable_reasons'].update(reference='invalid_or_missing_coordinates',before='invalid_or_missing_coordinates',after='invalid_or_missing_coordinates')
  sites.append(s)
 atomic(CACHE/'plan.json',{'reference_year':REF,'sites':sites,'jobs':jobs})
 return rows,sites,jobs
def seed(sites,jobs):
 # New land-coverage denominator: do not reuse the earlier geometric-coverage screening.
 return 0
def search(lon,lat,year):
 params={'collections':'landsat-c2-l2','bbox':f'{lon-.001},{lat-.001},{lon+.001},{lat+.001}','datetime':f'{year}-06-01T00:00:00Z/{year}-09-01T00:00:00Z','limit':'100'}
 j=get_json(BASE+'/search?'+urllib.parse.urlencode(params));items=j.get('features',[])
 while any(l.get('rel')=='next' for l in j.get('links',[])):
  nxt=next(l for l in j['links'] if l['rel']=='next')
  if nxt.get('method','GET')!='GET':raise RuntimeError('pagination_requires_unhandled_method')
  j=get_json(nxt['href']);items+=j.get('features',[])
 platforms=['landsat-8','landsat-9'] if year>=2013 else ['landsat-4','landsat-5','landsat-7']
 return sorted([f for f in items if f['properties'].get('platform') in platforms and f['properties'].get('landsat:correction')=='L2SP' and ('lwir11' in f['assets'] or 'lwir' in f['assets'])],key=lambda f:f['properties'].get('eo:cloud_cover',100))
def read(item,asset,lon,lat):
 if asset=='lwir11' and asset not in item['assets']:asset='lwir'
 url=item['assets'][asset]['href']+'?'+token()
 with rasterio.open(url) as ds:
  x,y=Transformer.from_crs(4326,ds.crs,always_xy=True).transform(lon,lat)
  window=from_bounds(x-8100,y-8100,x+8100,y+8100,ds.transform).round_offsets().round_lengths()
  a=ds.read(1,window=window,boundless=True,fill_value=0);tr=ds.window_transform(window)
  yy,xx=np.indices(a.shape);d=np.hypot(tr.c+(xx+.5)*tr.a-x,tr.f+(yy+.5)*tr.e-y)
  return a,d
def scene(item,lon,lat):
 qa,d=read(item,'qa_pixel',lon,lat);inner=d<1000;bg=(d>=6000)&(d<8000)
 land=(qa&128)==0;valid=(qa&255)==64
 def coverage(m):return float(np.sum(valid&m)/np.sum(land&m)) if np.sum(land&m) else 0.0
 def land_fraction(m):return float(np.mean(land[m]))
 def sufficient():return land_fraction(inner)>=.25 and land_fraction(bg)>=.1 and coverage(inner)>=.7 and coverage(bg)>=.7
 if not sufficient():return None
 unc,_=read(item,'qa',lon,lat);valid&=(unc>=0)&(unc*.01<=5)
 if not sufficient():return None
 dist,_=read(item,'cdist',lon,lat);valid&=dist*.01>=1
 if not sufficient():return None
 sat,_=read(item,'qa_radsat',lon,lat);valid&=sat==0
 if not sufficient():return None
 raw,_=read(item,'lwir11',lon,lat);valid&=raw>0
 if not sufficient():return None
 temp=np.where(valid,raw*.00341802+149-273.15,np.nan);background=float(np.nanmedian(temp[bg]));bands=[]
 for lo,hi in RINGS:
  m=(d>=lo)&(d<hi);cov=coverage(m);lf=land_fraction(m);v=float(np.nanmedian(temp[m])) if cov>=.7 and lf>=.25 else None
  bands.append({'inner_m':lo,'outer_m':hi,'coverage':cov,'land_fraction':lf,'lst_c':v,'background_c':background,'anomaly_c':v-background if v is not None else None})
 return {'scene_id':item['id'],'date':item['properties']['datetime'],'platform':item['properties']['platform'],'scene_cloud_pct':item['properties'].get('eo:cloud_cover'),'background_c':background,'inner_land_fraction':land_fraction(inner),'background_land_fraction':land_fraction(bg),'coverage_basis':'valid_land_pixels_over_nonwater_pixels','median_st_qa_k':float(np.median(unc[valid&inner])*.01),'fraction_inner_qa_le_2k':float(np.mean((unc[inner]*.01<=2)&valid[inner])),'bands':bands}
def run_job(job):
 path=PER/(job['key']+'.json')
 old=None
 if path.exists():
  old=json.loads(path.read_text(encoding='utf-8'))
  if old.get('status')!='access_error' and not (old.get('errors') and len(old.get('accepted',[]))<2):return old
 result={**job,'status':'access_error','accepted':[],'errors':[],'attempted':0,'available_scenes':0,'provenance':'direct_public_landsat_c2_l2'}
 if old:result.update(accepted=old.get('accepted',[]),retry_count=old.get('retry_count',0)+1)
 exhaustive=os.environ.get('HEAT_EXHAUSTIVE_ERRORS')=='1'
 excluded={'LC08_L2SP_015033_20220804_02_T1'} if exhaustive else set()
 if exhaustive:result.update(candidate_policy='exhaustive_same_window_after_access_failure',unreadable_qa_scenes_excluded=sorted(excluded))
 try:
  items=search(job['longitude'],job['latitude'],job['year']);result['available_scenes']=len(items)
  for item in (items if exhaustive else items[:8]):
   if item['id'] in excluded:continue
   if any(a['date'][:10]==item['properties']['datetime'][:10] for a in result['accepted']):continue
   result['attempted']+=1
   try:
    a=scene(item,job['longitude'],job['latitude'])
    if a:result['accepted'].append(a)
   except Exception as e:result['errors'].append({'scene_id':item['id'],'type':type(e).__name__,'http_code':getattr(e,'code',None),'step':'raster'})
   if len(result['accepted'])>=2:break
  result['status']='ok' if len(result['accepted'])>=2 else ('access_error' if result['errors'] else 'insufficient_scenes')
 except Exception as e:result['errors'].append({'type':type(e).__name__,'http_code':getattr(e,'code',None),'step':'catalog'})
 atomic(path,result)
 return result
def main():
 rows,sites,jobs=plan();seeded=seed(sites,jobs)
 todo=[j for k,j in jobs.items() if not (PER/(k+'.json')).exists()]
 if '--retry-errors' in sys.argv:
  todo=[]
  for k,j in jobs.items():
   if not (PER/(k+'.json')).exists():continue
   p=json.loads((PER/(k+'.json')).read_text(encoding='utf-8'))
   if p['status']=='access_error' or (p.get('errors') and len(p.get('accepted',[]))<2):todo.append(j)
 todo.sort(key=lambda j:(j['year']!=REF,j['year'],round(j['longitude'],1),j['latitude']))
 if '--reverse' in sys.argv:todo=list(reversed(todo))
 if '--probe' in sys.argv:todo=todo[:4]
 print(json.dumps({'all_records':len(rows),'valid_coordinates':sum(s['coordinate_valid'] for s in sites),'unique_coordinate_year_jobs':len(jobs),'reused_prior_periods':seeded,'remaining_jobs':len(todo),'reference_year':REF}),flush=True)
 counts=collections.Counter();start=time.time()
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.TIF,.tif',GDAL_HTTP_TIMEOUT='35',GDAL_HTTP_MAX_RETRY='1',CPL_VSIL_CURL_CACHE_SIZE=268435456,GDAL_CACHEMAX=268435456):
  with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.environ.get('HEAT_ALL_WORKERS','64'))) as pool:
   fs=[pool.submit(run_job,j) for j in todo]
   for n,f in enumerate(concurrent.futures.as_completed(fs),1):
    try:r=f.result();counts[r['status']]+=1
    except Exception as e:counts['worker_error']+=1;print(type(e).__name__,flush=True)
    if n%50==0 or n==len(todo):
     progress={'completed_this_run':n,'remaining_this_run':len(todo)-n,'statuses':dict(counts),'elapsed_minutes':round((time.time()-start)/60,1)};atomic(CACHE/('progress_reverse.json' if '--reverse' in sys.argv else 'progress_retry.json' if '--retry-errors' in sys.argv else 'progress.json'),progress);print(json.dumps(progress),flush=True)
 print('RUN_FINISHED',flush=True)
if __name__=='__main__':main()
