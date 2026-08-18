#!/usr/bin/env python3
import argparse
from pathlib import Path

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--samples", required=True)
    p.add_argument("--vcfs", nargs="+", required=True)
    p.add_argument("--output", required=True)
    a=p.parse_args()
    samples=[x for x in a.samples.split(",") if x]
    if len(samples) != len(a.vcfs):
        raise SystemExit(f"sample/VCF count mismatch: {len(samples)} != {len(a.vcfs)}")
    out=Path(a.output); out.parent.mkdir(parents=True, exist_ok=True)
    tmp=out.with_suffix(out.suffix + ".tmp")
    tmp.write_text("".join(f"{s}\t{Path(v).resolve()}\n" for s,v in zip(samples,a.vcfs)))
    tmp.replace(out)

if __name__ == "__main__": main()
