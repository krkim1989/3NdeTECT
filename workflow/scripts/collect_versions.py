#!/usr/bin/env python3
import argparse,platform,subprocess
from pathlib import Path
import pandas as pd

COMMANDS={"python":["python","--version"],"bwa":["bwa"],"samtools":["samtools","--version"],"bcftools":["bcftools","--version"],"gatk":["gatk","--version"],"nQuire":["nQuire"],"meryl":["meryl","--version"]}
def first_line(cmd):
    try:
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=30); text=(r.stdout+"\n"+r.stderr).strip(); return text.splitlines()[0] if text else f"exit={r.returncode}"
    except Exception as e: return f"unavailable: {type(e).__name__}"
def main():
    p=argparse.ArgumentParser(); p.add_argument("--output",required=True); a=p.parse_args(); rows=[{"software":k,"version":first_line(v)} for k,v in COMMANDS.items()]; rows.append({"software":"platform","version":platform.platform()})
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True); tmp=out.with_suffix(out.suffix+".tmp"); pd.DataFrame(rows).to_csv(tmp,sep="\t",index=False); tmp.replace(out)
if __name__=="__main__": main()
