"""Fresh pdflatex build; reject warnings, unresolved references and box defects."""
from pathlib import Path
import argparse,subprocess,shutil,re,sys,tempfile
p=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--build-dir',default=None,help='new or empty build directory; default: a new system temporary directory')
a=parser.parse_args();b=(p/a.build_dir).resolve() if a.build_dir else Path(tempfile.mkdtemp(prefix='exponential-ergodicity-build-')).resolve()
if b==p:parser.error('build directory must differ from the paper source directory')
if b.exists() and any(b.iterdir()):parser.error('build directory must be new or empty')
b.mkdir(parents=True,exist_ok=True)
def run(args,name):
 r=subprocess.run(args,cwd=p,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
 (b/name).write_text(r.stdout,encoding='utf-8',newline='\n')
 if r.returncode:print(r.stdout[-8000:]);sys.exit(r.returncode)
# Some LaTeX releases (for example the 2025-11-01 kernel) ship a first-aid patch for cleveref, the 'label' hook chunk
# 'firstaid/cleveref', which names every theorem-like environment sharing the theorem counter "Theorem". Removing it
# restores cleveref's own names (Proposition, Lemma, ...). Where the chunk is absent LaTeX only notes that fact.
# DoNotLoadEpstopdf keeps the graphics driver from loading epstopdf, which warns when shell escape is disabled.
bs=chr(92)
job=bs+'def'+bs+'DoNotLoadEpstopdf{}'+bs+'AddToHook{begindocument/end}{'+bs+'RemoveFromHook{label}[firstaid/cleveref]}'+bs+'input{main.tex}'
tex=['pdflatex','-no-shell-escape','-disable-installer','-recorder','-interaction=nonstopmode','-halt-on-error','-output-directory='+str(b),'-jobname=main',job]
for i in range(1,4):run(tex,'pdflatex-'+str(i)+'.stdout')
log=(b/'main.log').read_text(errors='replace')
bad=[l for l in log.splitlines() if re.search(r'Warning:|Overfull|Underfull|undefined references|undefined citations|ignored error',l,re.I) and "Cannot remove chunk 'firstaid/cleveref'" not in l]
if bad:print('FINAL BUILD ISSUES:\n'+'\n'.join(bad));sys.exit(2)
shutil.copy2(b/'main.pdf',p/'main.pdf')
print('Clean pdflatex build; zero warnings, unresolved references, overfull or underfull boxes.')
