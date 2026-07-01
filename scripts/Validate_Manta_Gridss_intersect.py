import pandas as pd
import re
from cyvcf2 import VCF, Writer
from collections import defaultdict
from sys import argv


def parse_Manta(manta_in, sid):
    linkdb = []
    infile = VCF(manta_in)
    event_dict = defaultdict(dict)
    for variant in infile:
        chrom = variant.CHROM
        pos = variant.POS
        mid = variant.ID
        event = mid.split(':', 1)[1].rsplit(':', 1)[0]
        svtype = variant.INFO.get("SVTYPE")
        if svtype in ["DEL","DUP","INV"]:
            end_pos = variant.INFO.get("END")
            #print(end_pos)
        else:
            end_pos = 0 
        
        # Store info grouped by EVENT
        count = len([k for k in event_dict[event].keys() if k.startswith("entry")])
        slot = f"entry{count+1}"
        event_dict[event][slot] = {
            "chrom": chrom,
            "pos": pos,
            "end_pos": end_pos,
            "svtype" : svtype,
            "variant": variant
        }

    for linkid, data in event_dict.items():
        row = {"sample": sid, "event": linkid}

        for i, (slot, record) in enumerate(data.items(), start=1):
            row[f"chrom_{i}"] = record["chrom"]
            row[f"pos_{i}"] = record["pos"]
            row[f"end_pos_{i}"] = record["end_pos"]
            row[f"type_{i}"] = record["svtype"]
            row[f"variant_{i}"] = record["variant"]
        linkdb.append(row)
     
    all_manta_variants = pd.DataFrame(linkdb)
    return all_manta_variants

def parse_Gridss(gridss_in, sampleID):
    gridss_event_dict = {}
    infile = VCF(gridss_in)
    
    for line in infile:
        gchrom = line.CHROM
        gpos = line.POS
        gevent = line.ID
        gsvtype = line.INFO.get("SVTYPE")
            
        event_entry = {
            "sample": sampleID,
            "event": gevent,
            "chrom": gchrom,
            "pos": gpos,
            "svtype": gsvtype,
            "variant": line  # keeping raw line if you need it later
        }
            
        # ensure list exists for sample
        gridss_event_dict.setdefault(sampleID, []).append(event_entry)
    all_events = [e for sample_events in gridss_event_dict.values() for e in sample_events]
    all_gridss_variants = pd.DataFrame(all_events)
    return all_gridss_variants

def define_overlap(manta_sub,gridss_sub,id):
    overlaps = []
    unmatched = []
    window = 50

    flank = 5000
    bed_lines = []
    for _, mrow in manta_sub.iterrows():
        matched_any = False
    
        for i in range(1, 5):  # up to 4 breakpoints
            m_chrom = mrow.get(f"chrom_{i}")
            m_pos = mrow.get(f"pos_{i}")
            
            if pd.isna(m_chrom) or pd.isna(m_pos):
                continue

            gridss_matches = gridss_sub[gridss_sub['chrom'] == m_chrom]
            if gridss_matches.empty:
                continue
            gridss_matches = gridss_matches.assign(
                distance=(gridss_matches['pos'] - m_pos).abs())    
            gridss_matches = gridss_matches[gridss_matches['distance'] <= window]

            if not gridss_matches.empty:
                closest = gridss_matches.loc[gridss_matches['distance'].idxmin()]
                overlaps.append({
                    "sample": id,
                    "sv_type": mrow[f'type_{i}'],
                    "end_pos" : mrow[f'end_pos_{i}'],
                    "manta_event": f"{mrow['event']}:{mrow[f'type_{i}']}",
                    "gridss_event": closest['event'],
                    "manta_chrom": m_chrom,
                    "manta_pos": m_pos,
                    "gridss_chrom": closest['chrom'],
                    "gridss_pos": closest['pos'],
                    "distance": closest['distance'],
                    "variant": mrow[f'variant_{i}']
                })
                matched_any = True
                
        if not matched_any:
            unmatched.append(mrow)
            v1 = mrow.get("variant_1")

            if v1 is not None:
                chrom = v1.CHROM
                pos = v1.POS
                start = max(0, pos - flank)
                end = pos + flank        # single bp interval
                bed_lines.append(f"{chrom}\t{start}\t{end}\t{mrow['event']}")

                
    if bed_lines:
        bed_file = f"{id}.gridss_manta_unmatched.bed"
        with open(bed_file,"w") as f:
            f.write("\n".join(bed_lines)+"\n")
        print(f"Saved {len(bed_lines)} breakpoints for sample {bed_file}")
            
    overlap = pd.DataFrame(overlaps)
    return overlap

def count_overlaps(all_df,sampleid):
    svtypes = ['DEL','INS','DUP','INV','BND']
    counts = (
        all_df.groupby(['sample','sv_type']).size().unstack(fill_value=0).reindex(columns=svtypes,fill_value=0)
    )
    counts.to_csv(sampleid + '.StructuralVariants.csv',sep="\t", index=True)
    return
    

sv_type_mapping = {
        'DEL': 'deletion',
        'INV': 'inversion',
        'DUP': 'tandem-duplication',
        'INS': 'insertion',
        'BND': 'translocation'
}

def write_to_bedpe(overlap_df,sample):
    bedpe_rows = []
    grouped = overlap_df.groupby("manta_event")  # or "event" if using simpler ID
    count = 0
    for event, group in grouped:
        svtype = group["sv_type"].iloc[0]
            
        if svtype in ["DEL","DUP","INV"]:
            chrom1 = group["manta_chrom"].iloc[0]
            start1 = group["manta_pos"].iloc[0] -1
            end1 = group["manta_pos"].iloc[0] 

            chrom2 = group["manta_chrom"].iloc[0]
            start2 = group["end_pos"].iloc[0] -1
            end2 = group["end_pos"].iloc[0]
            
            
        elif svtype == "BND":
            if len(group) == 2:
                chrom1 = group["manta_chrom"].iloc[0]
                start1 = group["manta_pos"].iloc[0] -1
                end1 = group["manta_pos"].iloc[0] 

                chrom2 = group["manta_chrom"].iloc[1]
                start2 = group["manta_pos"].iloc[1] -1
                end2 = group["manta_pos"].iloc[1]
            else:
                count=+1
                continue
                
        svclass = sv_type_mapping.get(svtype)
        bedpe_rows.append([
            sample, chrom1, start1, end1,
            chrom2, start2, end2, svclass, event
        ])

    bedpe_df = pd.DataFrame(
        bedpe_rows,
        columns=[
            "sample","chrom1", "start1", "end1",
            "chrom2", "start2", "end2", "svclass", "event"
        ]
    )
    print(f"for {sample}, {count} translocations are described with only 1 or more than 2 BND rows per manta event, they are excluded from the bedpe")
    bedpe_df.to_csv(f"{sample}.final.bedpe", sep="\t", index=False)
    return

def write_to_vcf(overlap_df,sample,input_vcf):
    inVCF = VCF(input_vcf)
    header = inVCF.raw_header
    nn = Writer(f"{sample}.structural.final.vcf", inVCF, "w")
    for variant in overlap_df['variant']:
        nn.write_record(variant)
    nn.close()
    return


manta_file=argv[2]
gridss_file=argv[3]
sample_id=argv[1]
manta_variants = parse_Manta(manta_file,sample_id)
gridss_variants = parse_Gridss(gridss_file,sample_id)
all_overlap = define_overlap(manta_variants,gridss_variants,sample_id)
write_to_bedpe(all_overlap,sample_id)
write_to_vcf(all_overlap,sample_id,manta_file)
all_overlap.to_csv(f"{sample_id}.all_structural_overlap.tsv", sep='\t', index=False)
overlap_counts = count_overlaps(all_overlap,sample_id)

