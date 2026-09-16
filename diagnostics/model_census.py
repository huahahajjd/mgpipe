"""Read MAT model dimensions without materializing model annotations."""
from pathlib import Path
import struct, zlib, json, os, time, concurrent.futures
BASE=Path('/mnt/nfs/wangchao/my_course/untreated_MS/replace_data/step4_modeling/predict_microbe_contributions')
OUT=BASE/'refactor_20260914'/'model_dimensions_v5.json'
class Stream:
 def __init__(self,iterator):self.it=iter(iterator);self.buf=b'';self.pos=0
 def take(self,n,keep=True):
  pieces=[];left=n
  while left:
   if not self.buf:self.buf=next(self.it)
   k=min(left,len(self.buf))
   if keep:pieces.append(self.buf[:k])
   self.buf=self.buf[k:];left-=k;self.pos+=k
  return b''.join(pieces) if keep else None
 def element(self):
  a,b=struct.unpack('<II',self.take(8));small=a>>16
  if small:return a&65535,small,struct.pack('<I',b)[:small]
  return a,b,None
 def value(self):
  t,n,inline=self.element()
  if inline is not None:return inline
  value=self.take(n);self.take((-n)%8,False);return value
 def skip(self,n):self.take(n,False)
def chunks(f,n):
 d=zlib.decompressobj();tail=b''
 while n or tail:
  if not tail:
   tail=f.read(min(n,1<<20));n-=len(tail)
   if not tail:break
  data=d.decompress(tail,1<<20);tail=d.unconsumed_tail
  if data:yield data
 data=d.flush()
 if data:yield data

def metadata(r,n):
 start=r.pos;flags=r.value();kind=struct.unpack('<I',flags[:4])[0]&255
 v=r.value();dims=struct.unpack('<'+'i'*(len(v)//4),v);name=r.value()
 return start,kind,dims,name

def top_model(r):
 typ,n,inline=r.element()
 if typ!=14:raise ValueError('Expected miMATRIX')
 start,kind,dims,name=metadata(r,n)
 if kind!=2:return None
 width=struct.unpack('<i',r.value())[0];fields=r.value();fields=[fields[i:i+width].split(b'\0')[0].decode() for i in range(0,len(fields),width)]
 if 'S' not in fields or 'rxns' not in fields:return None
 result={'couplingRows':0}
 for field in fields:
  typ,size,inline=r.element();assert typ==14 and inline is None
  st,cls,shape,_=metadata(r,size)
  if field=='S':result['reactions']=shape[1];result['metaboliteRows']=shape[0]
  if field=='C':result['couplingRows']=shape[0]
  if 'reactions' in result and ('C' not in fields or field=='C'):return result
  r.skip(size-(r.pos-st));r.skip((-size)%8)
 raise ValueError('S/C not found')
def mat5(f):
 with f.open('rb') as h:
  header=h.read(128);assert header[126:128]==b'IM'
  while True:
   raw=h.read(8)
   if not raw:break
   typ,n=struct.unpack('<II',raw);start=h.tell()
   if typ==15:r=Stream(chunks(h,n));result=top_model(r)
   elif typ==14:
    h.seek(start-8);r=Stream(iter(lambda:h.read(1<<20),b''));result=top_model(r)
   else:result=None
   if result:return result
   h.seek(start+n)
 raise ValueError('No model structure')
def inspect(f):
 st=f.stat();t=time.time()
 with f.open('rb') as h:header=h.read(128)
 assert b'MATLAB 7.3' not in header, 'Use census_hdf.m for v7.3 files'
 result=mat5(f)
 result.update(file=str(f),name=f.name,bytes=st.st_size,mtime=st.st_mtime,seconds=time.time()-t)
 return result
if __name__=='__main__':
 files=sorted((BASE/'run_20260507_204155'/'ordered_models').glob('*.mat'))
 records=json.loads(OUT.read_text()) if OUT.exists() else []
 known={r['name']:r for r in records};todo=[f for f in files if f.name not in known or known[f.name]['bytes']!=f.stat().st_size or known[f.name]['mtime']!=f.stat().st_mtime]
 def commit(result):
  known[result['name']]=result;tmp=OUT.with_suffix('.partial');tmp.write_text(json.dumps(list(known.values()),indent=2));os.replace(str(tmp),str(OUT));print(len(known),result['name'],result['reactions'],result['couplingRows'],flush=True)
 v5=[]
 for f in todo:
  with f.open('rb') as h:is73=b'MATLAB 7.3' in h.read(128)
  if is73:continue
  else:v5.append(f)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for result in pool.map(inspect,v5):commit(result)
 records=list(known.values());assert len(records)==352
 lowestR=min(r['reactions'] for r in records);lowestC=min(r['couplingRows'] for r in records)
 both=[r for r in records if r['reactions']==lowestR and r['couplingRows']==lowestC]
 print(json.dumps({'minReactions':lowestR,'minCouplingRows':lowestC,'jointMinimum':both,'smallestSum':min(records,key=lambda r:r['reactions']+r['couplingRows'])},indent=2),flush=True)
