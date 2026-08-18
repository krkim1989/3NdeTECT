from pathlib import Path
import yaml

path=Path(snakemake.output.yaml); path.parent.mkdir(parents=True,exist_ok=True)
tmp=path.with_suffix(path.suffix+".tmp"); tmp.write_text(yaml.safe_dump(dict(snakemake.config),sort_keys=False)); tmp.replace(path)
