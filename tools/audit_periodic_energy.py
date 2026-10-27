import argparse,csv,json,math
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description="Export saved periodic stability histories; does not run a solver.")
    p.add_argument("--input",type=Path,nargs="+",required=True)
    p.add_argument("--output",type=Path,default=Path("periodic-energy-export"))
    a=p.parse_args();files=[]
    for item in a.input:
        if item.is_file():files.append(item)
        elif item.is_dir():
            files.extend(item.rglob("stability_periodic*.csv"))
            files.extend(item.rglob("periodic_audit_*_history.csv"))
        else:raise ValueError(f"Missing input: {item}")
    files=sorted(set(x.resolve() for x in files))
    if not files:raise ValueError("No periodic time histories found; summary CSVs cannot reconstruct them")
    a.output.mkdir(parents=True,exist_ok=True); summaries=[];legacy=[]
    for path in files:
        with path.open(newline="") as f:rows=list(csv.DictReader(f))
        if not rows:raise ValueError(f"Empty history: {path}")
        modern="modified_energy_j" in rows[0]
        stepkey="completed_steps" if modern else "step"
        tkey="t_e_s" if modern else "t"
        energykey="raw_minus_energy_j" if modern else "diagnostic_energy"
        for row in rows:
            for key in [stepkey,tkey,energykey,"modified_energy_j" if modern else "interior_energy"]:
                if not math.isfinite(float(row[key])):raise ValueError(f"Nonfinite {key}: {path}")
            if float(row["source_active"])!=0:raise ValueError(f"Source-active row in periodic history: {path}")
        steps=[int(row[stepkey]) for row in rows]
        if any(y<=x for x,y in zip(steps,steps[1:])):raise ValueError(f"Non-increasing steps: {path}")
        energy=[float(row[energykey]) for row in rows];ref=energy[0]
        if ref<=0:raise ValueError(f"Nonpositive first-row reference: {path}")
        summary={"source_file":str(path),"rows":len(rows),"first_step":steps[0],"last_step":steps[-1],"sampled_min_over_first":min(energy)/ref,"sampled_max_over_first":max(energy)/ref,"final_over_first":energy[-1]/ref,"max_abs_modified_relative_drift":None}
        if modern:
            summary["max_abs_modified_relative_drift"]=max(abs(float(row["modified_relative_drift"])) for row in rows)
        else:
            for row in rows:legacy.append({"source_file":str(path),**row,"diagnostic_energy_over_first":float(row[energykey])/ref})
        summaries.append(summary)
        target=a.output/(path.stem+"_export.csv")
        if target.exists():raise ValueError(f"Refusing to overwrite: {target}")
        with target.open("w",newline="") as f:
            fields=list(rows[0]);w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    if legacy:
        with (a.output/"periodic_legacy_normalized.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(legacy[0]));w.writeheader();w.writerows(legacy)
    with (a.output/"periodic_energy_summary.json").open("w") as f:json.dump(summaries,f,indent=2,allow_nan=False)
    print(json.dumps(summaries,indent=2,allow_nan=False))
if __name__=="__main__":main()