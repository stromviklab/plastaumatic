#!/usr/bin/env python3

import sys 
import argparse
import pathlib
import csv
from datetime import date, datetime
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline

## parse the arguments 
parser=argparse.ArgumentParser()
parser.add_argument('-f','--fasta', metavar='\b', type=pathlib.Path, required=True, help="input fasta file")
parser.add_argument('-g','--genbank', metavar='\b', type=pathlib.Path, required=True, help="reference genbank file")
parser.add_argument('-o','--out', metavar='\b', type=pathlib.Path, default="out_target", help="output directory")
parser.add_argument('-p','--prefix', metavar='\b', default="target", help="prefix for the output files")
args=parser.parse_args() 

## check if input files exist
if not (args.fasta.absolute().is_file() and args.genbank.absolute().is_file()):
    print("ERROR:","Input fasta/genbank files does not exist", sep='\t')
    sys.exit()

## create output directory if does not exist 
args.out.absolute().mkdir(parents=True, exist_ok=True)

## write reference features to a file
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Extracting referene features", sep='\t')
ref_seqs=open(args.out.absolute()/"ref_features.fa",'w')
my_fa=SeqIO.read(args.fasta,"fasta")
my_gb=SeqIO.read(args.genbank,"genbank")
ftr_types="gene CDS intron exon rRNA tRNA"
smallest_ftr=[]
ftr_ids=[]
for feature in my_gb.features[1:]:
    if any(ftr_type == feature.type for ftr_type in ftr_types.split()):
        for part in feature.location.parts:
            if "gene" in feature.qualifiers:
                print(">",list(feature.qualifiers["gene"])[0],"_",feature.type,"_",part.nofuzzy_start,"_",part.nofuzzy_end,sep="", file=ref_seqs)
                print(part.extract(my_gb.seq), file=ref_seqs)
                smallest_ftr.append(part.nofuzzy_end-part.nofuzzy_start)
                ftr_ids.append(list(feature.qualifiers["gene"])[0]+"_"+feature.type+"_"+str(part.nofuzzy_start)+"_"+str(part.nofuzzy_end))
            else:
                print(">id_",feature.type,"_",part.nofuzzy_start,"_",part.nofuzzy_end,sep="", file=ref_seqs)
                print(part.extract(my_gb.seq), file=ref_seqs)
                smallest_ftr.append(part.nofuzzy_end-part.nofuzzy_start)
                ftr_ids.append("id_"+feature.type+"_"+str(part.nofuzzy_start)+"_"+str(part.nofuzzy_end))
ref_seqs.close()
ftr_ids=set(ftr_ids)

## blast alignment and filtering function
def run_blast(p_ident):
    if min(smallest_ftr) > 10:
        blast_cmd = NcbiblastnCommandline(task="blastn", query=args.out.absolute()/"ref_features.fa", subject=args.fasta.absolute(), \
            perc_identity=p_ident, max_hsps=5, max_target_seqs=5, evalue=100, out=args.out.absolute()/"blast.anno.out", outfmt="6 qseqid qstart qend sseqid sstart send length pident sstrand qcovhsp")
    else:
        blast_cmd = NcbiblastnCommandline(task="blastn", query=args.out.absolute()/"ref_features.fa", subject=args.fasta.absolute(), \
            perc_identity=p_ident, max_hsps=50, max_target_seqs=5, evalue=1000, word_size=min(smallest_ftr), out=args.out.absolute()/"blast.anno.out", outfmt="6 qseqid qstart qend sseqid sstart send length pident sstrand qcovhsp")
    blast_cmd()

    ## Parse blast output 
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Parsing the blast output", sep='\t')
    pd.options.mode.chained_assignment = None
    df=pd.read_table(args.out.absolute()/"blast.anno.out",sep='\t',names=["ref","r_start","r_end","qry","q_start","q_end","evalue","ident","strand","cov"])

    ### filter alignments and separate into duplicated and non-duplicated alignemnts 
    df2=df.loc[df['cov']>99]
    df2=df2[['ref','r_start','r_end','q_start','q_end','strand']]
    df_plus=df2.loc[df2['strand'] == "plus"]
    df_minus=df2.loc[df2['strand'] == "minus",['ref','r_start','r_end','q_end','q_start','strand']]
    df_minus.rename(columns={"q_start":"q_end","q_end":"q_start"},inplace=True)
    df2=pd.concat([df_plus,df_minus])
    df2=df2.loc[~df2.duplicated(['ref','q_start','q_end'], keep="last")]

    df_uniq=df2.loc[~df2.duplicated('ref',keep=False)].copy()
    df_dup=df2.loc[df2.duplicated('ref',keep=False)].copy()
    cols=['ref','type','r_start','r_end','q_start','q_end','strand']
    new1 = df_uniq["ref"].str.split("_",expand=True) 
    new1[2]=pd.to_numeric(new1[2])
    new1[3]=pd.to_numeric(new1[3])
    df_uniq=new1.join(df_uniq) 
    df_uniq.drop(['ref','r_start','r_end'],axis=1, inplace=True)
    df_uniq.columns=cols
    new2 = df_dup["ref"].str.split("_",expand=True)
    new2[2]=pd.to_numeric(new2[2])
    new2[3]=pd.to_numeric(new2[3])
    df_dup=new2.join(df_dup)
    df_dup.drop(['ref','r_start','r_end'],axis=1, inplace=True)
    df_dup.columns=cols

    ### get all features from non-duplicated alignments 
    gene_uniq=df_uniq.loc[df_uniq['type']=='gene']
    mergedStuff = pd.merge(gene_uniq, df_dup, on=['ref'], how='inner')
    map=mergedStuff[(mergedStuff.r_start_y >= mergedStuff.r_start_x) & (mergedStuff.r_end_y <= mergedStuff.r_end_x) & (mergedStuff.q_start_y >= mergedStuff.q_start_x) & (mergedStuff.q_end_y <= mergedStuff.q_end_x)]
    map.drop(['type_x','r_start_x','r_end_x','q_start_x','q_end_x','strand_x'],axis=1,inplace=True)
    map.columns=cols 
    gene_map1=pd.concat([df_uniq,map])

    ### drop features from duplicated if they are already written 
    drop_df=pd.merge(gene_map1, df_dup, on=['ref','type','r_start','r_end'], how='inner',suffixes=("_x",None))
    drop_df.drop(['q_start_x','q_end_x','strand_x'],axis=1,inplace=True)
    df_dup2=pd.concat([df_dup,drop_df]).drop_duplicates(keep=False)

    ### first, get all gene features from duplicated alignments 
    if len(df_dup2) > 0:
        df_dup2_gene=df_dup2.loc[df_dup2['type']=='gene','ref'].drop_duplicates(keep="first")
        dup2_gene_final=[]
        for i in df_dup2_gene:
            tmp=df_dup2.loc[(df_dup2['ref']==i) & (df_dup2['type']=='gene')]
            tmp1=tmp.loc[tmp.duplicated(['ref','type','r_start','r_end'],keep="first")].copy()
            tmp2=tmp.loc[~tmp.duplicated(['ref','type','r_start','r_end'],keep="first")].copy()
            tmp_merged=pd.merge(tmp1,tmp2,how="inner",on=['ref','type','r_start','r_end'])
            if len(tmp_merged) >= 2:
                for ind, row in tmp_merged.iterrows():
                    row2=row.tolist()
                    if ind == 0:
                        dup2_gene_final.append(row2[0:7])
                    # -10 is added for some features where the qry positions are not exactly same
                    elif any(lst[4] >= row2[4]-10 and lst[5] <= row2[5]+10 for lst in dup2_gene_final):  
                        dup2_gene_final.append(row2[0:4]+row2[7:10])
                    else:
                        dup2_gene_final.append(row2[0:7])
            else:
                tmp_nocds=pd.merge(tmp_merged, gene_map1.loc[gene_map1['type']=='gene'], on=['ref'], how='inner')
                row2=tmp_nocds.loc[0].tolist()
                if (~(tmp_nocds.q_start_x >= tmp_nocds.q_start) & (tmp_nocds.q_end_x <= tmp_nocds.q_end)).bool():
                    dup2_gene_final.append(row2[0:7])
                else:
                    dup2_gene_final.append(row2[0:4]+row2[7:10]) 
        dup2_gene_final_df=pd.DataFrame(dup2_gene_final)
        dup2_gene_final_df.columns=cols   

        ### get all non-gene features from duplicated alignments 
        nongene_dup2=df_dup2.loc[df_dup2['type']!='gene']
        mergedStuff2 = pd.merge(dup2_gene_final_df, nongene_dup2, on=['ref'], how='inner',suffixes=("_x",None))
        map2=mergedStuff2[(mergedStuff2.r_start >= mergedStuff2.r_start_x) & (mergedStuff2.r_end <= mergedStuff2.r_end_x) & (mergedStuff2.q_start >= mergedStuff2.q_start_x) & (mergedStuff2.q_end <= mergedStuff2.q_end_x)]
        map2.drop(['type_x','r_start_x','r_end_x','q_start_x','q_end_x','strand_x'],axis=1,inplace=True)
        map2.columns=cols 
        gene_map2=pd.concat([dup2_gene_final_df,map2])

        ### concatenate all the mappings
        final_gene_map=pd.concat([gene_map1,gene_map2])
    else:
        final_gene_map=gene_map1
    final_gene_map['q_start'] -= 1
    return final_gene_map
    

## round 1 - percent identity of 95 
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Running round-1 blastn searches", sep='\t')
final_gene_map=run_blast(95)
if len(ftr_ids) == len(final_gene_map):
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Found all the features from the reference", sep='\t')    
else:
    ## round 2 - percent identity of 90 
    print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Running round-2 blastn searches", sep='\t')
    final_gene_map=run_blast(90)
    if len(ftr_ids) == len(final_gene_map):
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Found all the features from the reference", sep='\t')    
    else:
        ## round 3 - percent identity of 85
        print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Running round-3 blastn searches", sep='\t')
        final_gene_map=run_blast(85)
        if len(ftr_ids) == len(final_gene_map):
            print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Found all the features from the reference", sep='\t') 
        else:
            print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"WARNING: some features are missing, use a different reference genbank file or annotate manually", sep='\t') 
final_gene_map.to_csv(args.out.absolute()/"qry-ref-gene.map",sep='\t',header=False, index=False)

## read the reference to target mappings
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Generating the target genebank file", sep='\t')
map_file=[]
with open(args.out.absolute()/"qry-ref-gene.map",newline='\n') as f:
    for row in csv.reader(f, delimiter='\t'):
        map_file.append(row)

## modify the referene features with the target  
new_features=[]
for record in SeqIO.parse(args.genbank,"genbank"):
    record.seq=my_fa.seq
    record.id=my_fa.id
    record.name=my_fa.name
    record.description=my_fa.id+" chloroplast, complete genome"
    record.annotations.clear()
    new_date=date.today().strftime("%d-%^b-%Y")
    new_record_anno={
        "molecule_type":"DNA",
        "topology":"circular",
        "data_file_division":"PLN",
        "date":new_date,
        "accessions":my_fa.id,
        "organism":my_fa.id,
        "sequence_version":1,
    }
    record.annotations.update(new_record_anno)
    record.dbxrefs.clear()    

    for feature in record.features[0:1]:
        feature.location=FeatureLocation(0,len(my_fa),1)
        feature.qualifiers.clear()
        new_qualifiers={
            "organelle":"plastid:chloroplast",
            "mol_type":"genomic DNA"
        }
        feature.qualifiers.update(new_qualifiers)
        new_features.append(feature)
    for feature in record.features[1:]:
        feature.qualifiers.pop('protein_id',None)
        feature.qualifiers.pop('db_xref',None)
        feature.qualifiers.pop('locus_tag',None)
        if len(feature.location.parts)==1:
            for my_list in map_file:
                if "gene" in feature.qualifiers:
                    if all(ftr in my_list[0:4] for ftr in (feature.qualifiers["gene"][0], feature.type, str(feature.location.nofuzzy_start), str(feature.location.nofuzzy_end))):
                        new_strand=1 if my_list[6]=="plus" else -1
                        new_loc=FeatureLocation(int(my_list[4]),int(my_list[5]),new_strand)
                        feature.location=new_loc
                        new_features.append(feature)
                        if "translation" in feature.qualifiers:
                            tr_seq=feature.location.extract(my_fa.seq)
                            tr_seq2=tr_seq.translate(table=11)
                            new_trnsl={"translation":tr_seq2}
                            feature.qualifiers.update(new_trnsl)
                            break
                        else:
                            break
                else:
                    if all(ftr in my_list[0:4] for ftr in ("id", feature.type, str(feature.location.nofuzzy_start), str(feature.location.nofuzzy_end))):
                        new_strand=1 if my_list[6]=="plus" else -1
                        new_loc=FeatureLocation(int(my_list[4]),int(my_list[5]),new_strand)
                        feature.location=new_loc
                        new_features.append(feature)
                        if "translation" in feature.qualifiers:
                            tr_seq=feature.location.extract(my_fa.seq)
                            tr_seq2=tr_seq.translate(table=11)
                            new_trnsl={"translation":tr_seq2}
                            feature.qualifiers.update(new_trnsl)
                            break
                        else:
                            break                        

        elif len(feature.location.parts)>1:
            i_count=0
            for i in range(0,len(feature.location.parts)):
                part=feature.location.parts[i]            
                for my_list in map_file:
                    if "gene" in feature.qualifiers:                        
                        if all(ftr in my_list[0:4] for ftr in (feature.qualifiers["gene"][0], feature.type, str(part.nofuzzy_start), str(part.nofuzzy_end))):
                            new_strand=1 if my_list[6]=="plus" else -1
                            new_loc=FeatureLocation(int(my_list[4]),int(my_list[5]),new_strand)
                            feature.location.parts[i]=new_loc
                            i_count=i_count+1
                            if i==len(feature.location.parts)-1 and i_count==len(feature.location.parts):
                                new_features.append(feature)
                                if "translation" in feature.qualifiers:
                                    tr_seq=feature.location.extract(my_fa.seq)
                                    tr_seq2=tr_seq.translate(table=11)
                                    new_trnsl={"translation":tr_seq2}
                                    feature.qualifiers.update(new_trnsl)
                                    break
                                else:
                                    break                                    
                    else:
                        if all(ftr in my_list[0:4] for ftr in ("id", feature.type, str(part.nofuzzy_start), str(part.nofuzzy_end))):
                            new_strand=1 if my_list[6]=="plus" else -1
                            new_loc=FeatureLocation(int(my_list[4]),int(my_list[5]),new_strand)
                            feature.location.parts[i]=new_loc
                            i_count=i_count+1
                            if i==len(feature.location.parts)-1 and i_count==len(feature.location.parts):
                                new_features.append(feature)
                                if "translation" in feature.qualifiers:
                                    tr_seq=feature.location.extract(my_fa.seq)
                                    tr_seq2=tr_seq.translate(table=11)
                                    new_trnsl={"translation":tr_seq2}
                                    feature.qualifiers.update(new_trnsl)   
                                    break
                                else:
                                    break                                    
                                                  
    record.features=new_features
    with open(args.out.absolute()/(args.prefix+'.gb'), "w") as new_gb:
        SeqIO.write(record, new_gb, "genbank")

print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"Finished", sep='\t')

