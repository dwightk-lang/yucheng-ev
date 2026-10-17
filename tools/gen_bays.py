"""產生 bays.json:以固定點(南港路三段47巷55號)為中心,把台北市路邊停車格位的靜態座標(WGS84 中心點)
連同「路段代碼 rdcode + 格位編號 pkid」存成精簡 JSON,供 roadside.html 與停管處即時檔 TCMSV_roadquery.xml 的
<roadSegID>/<psId> 對應(即時檔本身座標全為 0,必須靠這份靜態對照)。

資料來源(shapefile, TWD97 TM2 / EPSG:3826, 約每月更新):
  data.taipei 資料集「臺北市路邊停車格位」(id 5a911ea5-1694-4301-808e-e1780d971611)
  https://data.taipei/api/dataset/5a911ea5-1694-4301-808e-e1780d971611/resource/7e2f32a0-9201-4666-a0ed-11034e6c2b66/download
  (zip 約 56MB,解開 ~525MB;內含 .shp/.shx/.dbf/.prj/.cpg)

用法:
  python3 -I tools/gen_bays.py <解壓目錄> bays.json [itp_near.json] [xml_cells_near.json]
  後兩個為選用的驗證輸入(itaipeiparking GetParks SelectType=5 的格位中心點 / XML 鄰近路段的 (roadSegID,psId) 清單)。
  輸出只保留 rdcode 與 pkid 皆非空的格位(可對應即時狀態者);2026-10-17 產出 4010 格、158 路段、約 178KB。
pktype 代碼(由 pktype1 文字反推):01 小型車、03 身障(汽車)、07/70 貨車裝卸、10 時段性、02 機車、22 汽機車彈性共用。
"""
import struct, math, json, sys, glob, collections, os
SHPDIR, OUT = sys.argv[1], sys.argv[2]
LAT0, LON0 = 25.051989, 121.591180      # 南港路三段47巷55號
RADIUS = 1600.0                          # m, preload radius
# --- TWD97 TM2 (EPSG:3826) <-> WGS84 (GRS80, k0=0.9999, lon0=121, FE=250000) ---
a=6378137.0; f=1/298.257222101; k0=0.9999; lon0=math.radians(121.0); FE=250000.0
e2=f*(2-f); ep2=e2/(1-e2)
def fwd(lat,lon):
    phi=math.radians(lat); lam=math.radians(lon)-lon0
    N=a/math.sqrt(1-e2*math.sin(phi)**2); T=math.tan(phi)**2; C=ep2*math.cos(phi)**2; A=lam*math.cos(phi)
    M=a*((1-e2/4-3*e2**2/64-5*e2**3/256)*phi-(3*e2/8+3*e2**2/32+45*e2**3/1024)*math.sin(2*phi)+(15*e2**2/256+45*e2**3/1024)*math.sin(4*phi)-(35*e2**3/3072)*math.sin(6*phi))
    x=FE+k0*N*(A+(1-T+C)*A**3/6+(5-18*T+T*T+72*C-58*ep2)*A**5/120)
    y=k0*(M+N*math.tan(phi)*(A*A/2+(5-T+9*C+4*C*C)*A**4/24+(61-58*T+T*T+600*C-330*ep2)*A**6/720))
    return x,y
def inv(x,y):
    x-=FE; M=y/k0; mu=M/(a*(1-e2/4-3*e2**2/64-5*e2**3/256)); e1=(1-math.sqrt(1-e2))/(1+math.sqrt(1-e2))
    phi1=mu+(3*e1/2-27*e1**3/32)*math.sin(2*mu)+(21*e1**2/16-55*e1**4/32)*math.sin(4*mu)+(151*e1**3/96)*math.sin(6*mu)+(1097*e1**4/512)*math.sin(8*mu)
    N1=a/math.sqrt(1-e2*math.sin(phi1)**2); T1=math.tan(phi1)**2; C1=ep2*math.cos(phi1)**2
    R1=a*(1-e2)/(1-e2*math.sin(phi1)**2)**1.5; D=x/(N1*k0)
    lat=phi1-(N1*math.tan(phi1)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*ep2)*D**4/24+(61+90*T1+298*C1+45*T1*T1-252*ep2-3*C1*C1)*D**6/720)
    lon=lon0+(D-(1+2*T1+C1)*D**3/6+(5-2*C1+28*T1-3*C1*C1+8*ep2+24*T1*T1)*D**5/120)/math.cos(phi1)
    return math.degrees(lat), math.degrees(lon)
X0,Y0=fwd(LAT0,LON0)
print('center TWD97', round(X0,2), round(Y0,2), 'roundtrip', inv(X0,Y0))
# --- DBF ---
dbf=glob.glob(os.path.join(SHPDIR,'*.dbf'))[0]; shp=glob.glob(os.path.join(SHPDIR,'*.shp'))[0]; shx=glob.glob(os.path.join(SHPDIR,'*.shx'))[0]
fd=open(dbf,'rb'); hdr=fd.read(32); nrec,hlen,rlen=struct.unpack('<xxxxIHH20x',hdr)
fields=[]; 
while True:
    b=fd.read(32)
    if b[0]==0x0D: break
    fields.append((b[:11].split(b'\0')[0].decode('ascii','ignore'),chr(b[11]),b[16]))
fpos={}; p=1
for n,t,l in fields: fpos[n]=(p,l); p+=l
def dbf_rec(i):
    fd.seek(hlen+i*rlen); rec=fd.read(rlen)
    g=lambda n:rec[fpos[n][0]:fpos[n][0]+fpos[n][1]].decode('utf-8','ignore').strip()
    return {n:g(n) for n in ['rdcode','pkid','pktype','pktype1','pkkind','roadname','pkroad','zone','area_name','pknum']}
# --- SHX/SHP ---
sx=open(shx,'rb'); sx.seek(100); idx=sx.read(); n_idx=len(idx)//8
sp=open(shp,'rb'); sp.seek(32); shptype=struct.unpack('<i',sp.read(4))[0]; print('shape type',shptype,'records(shx)',n_idx,'dbf',nrec)
bays=[]; types=collections.Counter(); kinds=collections.Counter(); bbox_hit=0
for i in range(n_idx):
    off,ln=struct.unpack('>ii',idx[i*8:i*8+8]); sp.seek(off*2+8)
    st=struct.unpack('<i',sp.read(4))[0]
    if st not in (5,15,25): continue
    xmin,ymin,xmax,ymax=struct.unpack('<4d',sp.read(32))
    if xmax<X0-RADIUS or xmin>X0+RADIUS or ymax<Y0-RADIUS or ymin>Y0+RADIUS: continue
    bbox_hit+=1
    nparts,npts=struct.unpack('<2i',sp.read(8)); parts=struct.unpack('<%di'%nparts,sp.read(4*nparts))
    pts=struct.unpack('<%dd'%(2*npts),sp.read(16*npts))
    xs=pts[0::2]; ys=pts[1::2]
    # centroid of first ring (drop closing point)
    end=parts[1] if nparts>1 else npts
    rx=xs[:end]; ry=ys[:end]
    if len(rx)>1 and rx[0]==rx[-1] and ry[0]==ry[-1]: rx=rx[:-1]; ry=ry[:-1]
    cx=sum(rx)/len(rx); cy=sum(ry)/len(ry)
    d=math.hypot(cx-X0,cy-Y0)
    if d>RADIUS: continue
    lat,lon=inv(cx,cy); r=dbf_rec(i)
    types[(r['pktype'],r['pktype1'])]+=1; kinds[r['pkkind']]+=1
    bays.append({'rd':r['rdcode'],'ps':r['pkid'],'lat':round(lat,6),'lon':round(lon,6),'name':r['roadname'] or r['pkroad'],'type':r['pktype'],'kind':r['pkkind'],'d':round(d)})
print('bbox hits',bbox_hit,'bays within radius',len(bays))
print('pktype counts',dict(types.most_common(12)))
print('pkkind counts',dict(kinds.most_common(8)))
print('empty rdcode/pkid',sum(1 for b in bays if not b['rd'] or not b['ps']))
# --- validate against itaipeiparking wkt centroids (itp_near.json: name/parkId) ---
import re
near=json.load(open(sys.argv[3])) if len(sys.argv)>3 else []
def norm(s):
    s=re.sub(r'^\S*\(\d+\)\s*','',s or '').strip()
    for a_,b_ in [('一段','1段'),('二段','2段'),('三段','3段'),('四段','4段'),('五段','5段'),('六段','6段')]: s=s.replace(a_,b_)
    return s
byk=collections.defaultdict(list)
for b in bays: byk[(norm(b['name']),b['ps'])].append(b)
errs=[]; nm=0
for it in near:
    c=byk.get((norm(it['name']),str(it['parkId'])))
    if not c: nm+=1; continue
    b=min(c,key=lambda b:math.hypot((b['lat']-it['lat'])*110540,(b['lon']-it['lon'])*111320*math.cos(math.radians(LAT0))))
    errs.append(math.hypot((b['lat']-it['lat'])*110540,(b['lon']-it['lon'])*111320*math.cos(math.radians(LAT0))))
if errs:
    errs.sort(); print(f'validation vs itaipeiparking wkt: matched {len(errs)}/{len(near)} (no-name-match {nm}); median err {errs[len(errs)//2]:.1f} m, p90 {errs[int(len(errs)*.9)]:.1f} m, max {errs[-1]:.1f} m')
# --- join coverage with XML cells ---
xmlcells=json.load(open(sys.argv[4])) if len(sys.argv)>4 else None
if xmlcells:
    keys=set((b['rd'],b['ps']) for b in bays); hit=sum(1 for k in xmlcells if tuple(k) in keys)
    print(f'XML cells of nearby segments found in shapefile: {hit}/{len(xmlcells)}')
    miss=[k for k in xmlcells if tuple(k) not in keys][:10]; print('sample missing',miss)
# names per rdcode
segs={}
for b in bays: segs.setdefault(b['rd'],collections.Counter())[b['name']]+=1
bays=[b for b in bays if b['rd'] and b['ps']]   # 只留可與即時檔對應的格位
segs={k:v for k,v in segs.items() if any(b['rd']==k for b in bays)}
out={'center':[LAT0,LON0],'radius':RADIUS,'source':'data.taipei 臺北市路邊停車格位 shapefile '+os.path.basename(shp),'bays':[[b['rd'],b['ps'],b['lat'],b['lon'],b['type']] for b in sorted(bays,key=lambda b:b['d'])],'segs':{k:v.most_common(1)[0][0] for k,v in segs.items()}}
json.dump(out,open(OUT,'w',encoding='utf-8'),ensure_ascii=False,separators=(',',':'))
print('wrote',OUT,os.path.getsize(OUT),'bytes; segments',len(segs))
print('sample',out['bays'][:3])
