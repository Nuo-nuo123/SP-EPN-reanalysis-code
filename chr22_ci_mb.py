# chr22 OR CI + per-Mb normalization sensitivity
import csv
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None

# DMP chromosome counts
dmp = list(csv.DictReader(open(r'D:\claude-study\Spinal Ependymoma\analysis\bulk\dmp_significant_annotated.csv', encoding='utf-8')))
from collections import Counter
chr_dmp = Counter(r['CHR'] for r in dmp)
print('DMPs:', len(dmp), '| chr22 DMPs:', chr_dmp.get('22', chr_dmp.get('chr22', 0)))

# background: deposited universe chromosome counts
mat = []
with open(r'D:\claude-study\Spinal Ependymoma\analysis\bulk\beta_72samples.csv', encoding='utf-8') as fh:
    rd = csv.reader(fh); next(rd)
    for r in rd:
        mat.append(r[0])
print('matrix probes:', len(mat))
# manifest CHR per probe
man = {}
with open(r'D:\claude-study\Spinal Ependymoma\analysis\bulk\epic_manifest_full.csv', encoding='utf-8') as fh:
    rd = csv.DictReader(fh)
    for r in rd:
        man[r['IlmnID']] = r['CHR']
chr_bg = Counter(man.get(p) for p in mat)
print('bg chr22:', chr_bg.get('22'), '| total with chr:', sum(1 for v in chr_bg.values() if v))

# OR + CI via normal approx (log OR)
import math
a = chr_dmp.get('22', 0)          # DMPs on chr22
b = len(dmp) - a                  # DMPs elsewhere
c = chr_bg.get('22', 0) - a       # non-DMP background chr22
d = (len(mat) - chr_bg.get('22', 0)) - b
print('2x2:', a, b, c, d)
logor = math.log(a*d/(b*c)) if a*b*c*d else None
se = math.sqrt(1/a + 1/b + 1/c + 1/d)
lo = math.exp(logor - 1.96*se); hi = math.exp(logor + 1.96*se)
print('OR = %.2f (95%% CI %.2f-%.2f)' % (math.exp(logor), lo, hi))

# per-Mb normalization: need chromosome lengths (Mb)
chr_len = {'1':249.25,'2':243.20,'3':198.02,'4':191.15,'5':180.92,'6':171.12,'7':159.14,'8':146.36,
           '9':141.21,'10':135.53,'11':135.01,'12':133.85,'13':115.17,'14':107.35,'15':102.53,
           '16':90.35,'17':83.26,'18':80.37,'19':58.62,'20':64.44,'21':46.71,'22':50.82}
dmp_mb = Counter()
bg_mb = Counter()
for p in mat:
    ch = man.get(p)
    if ch in chr_len:
        bg_mb[ch] += 1
for r in dmp:
    ch = r['CHR']
    if ch in chr_len:
        dmp_mb[ch] += 1
# DMPs per Mb per chromosome
rows = []
for ch, ln in chr_len.items():
    if bg_mb.get(ch, 0) > 0:
        rows.append((ch, dmp_mb.get(ch, 0)/ln, bg_mb[ch]/ln, (dmp_mb.get(ch, 0)/ln)/(bg_mb[ch]/ln)))
rows.sort(key=lambda x: -x[0])
print('\nDMPs per Mb (top 5 + chr22):')
for ch, dm, bm, ratio in rows[:5] + [x for x in rows if x[0] == '22']:
    print('  chr%s: DMP/Mb %.2f | probes/Mb %.0f | ratio %.2f' % (ch, dm, bm, ratio))
ch22 = [x for x in rows if x[0] == '22'][0]
med_ratio = sorted(x[3] for x in rows)[len(rows)//2]
print('chr22 ratio %.2f vs median chr ratio %.2f' % (ch22[3], med_ratio))
