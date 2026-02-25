from __future__ import annotations
import numpy as np
import pandas as pd
from itertools import combinations

def _posterior_mdr_prob(row: pd.Series, threshold:int=3, missing_prior: dict|None=None) -> float:
    vals=row.to_dict(); known_pos=sum(1 for v in vals.values() if pd.notna(v) and int(v)==1); missing=[k for k,v in vals.items() if pd.isna(v)]
    if known_pos>=threshold: return 1.0
    if known_pos+len(missing)<threshold: return 0.0
    p=0.0
    for mask in range(1<<len(missing)):
        kk=known_pos; pr=1.0
        for i,c in enumerate(missing):
            pc=float((missing_prior or {}).get(c,0.5)); bit=(mask>>i)&1
            if bit: kk+=1; pr*=pc
            else: pr*=(1-pc)
        if kk>=threshold: p+=pr
    return float(min(max(p,0.0),1.0))

def next_best_test_evpi(class_df: pd.DataFrame, mdr_threshold: int=3, test_costs: dict|None=None, fn_cost: float=5.0, fp_cost: float=1.0) -> pd.DataFrame:
    if class_df.empty: return pd.DataFrame()
    prev=class_df.mean(skipna=True).fillna(0.5).to_dict(); costs={c:1.0 for c in class_df.columns}
    if test_costs:
        for k,v in test_costs.items():
            if k in costs:
                try: costs[k]=max(float(v),1e-9)
                except Exception: pass
    rows=[]
    for idx,row in class_df.iterrows():
        miss=[c for c in class_df.columns if pd.isna(row[c])]
        if not miss: continue
        p_mdr=_posterior_mdr_prob(row, threshold=mdr_threshold, missing_prior=prev)
        risk0=min(p_mdr*fn_cost, (1-p_mdr)*fp_cost)
        for c in miss:
            pc=float(prev.get(c,0.5))
            r1=row.copy(); r1[c]=1
            r0=row.copy(); r0[c]=0
            p1=_posterior_mdr_prob(r1, threshold=mdr_threshold, missing_prior=prev)
            p0=_posterior_mdr_prob(r0, threshold=mdr_threshold, missing_prior=prev)
            risk_after = pc*min(p1*fn_cost, (1-p1)*fp_cost) + (1-pc)*min(p0*fn_cost, (1-p0)*fp_cost)
            evpi=max(risk0-risk_after,0.0)
            rows.append({'Strain_ID':idx,'Candidate_Test_Class':c,'Current_pMDR':p_mdr,'Current_BayesRisk':risk0,'Expected_BayesRisk_AfterTest':risk_after,'EVPI_RiskReduction':evpi,'Test_Cost':costs[c],'EVPI_per_Cost':evpi/costs[c],'Class_Prevalence':pc,'Current_Uncertain':int(0<p_mdr<1),'N_tested':int(row.notna().sum())})
    out=pd.DataFrame(rows)
    return out.sort_values(['Strain_ID','EVPI_per_Cost','EVPI_RiskReduction'], ascending=[True,False,False]) if not out.empty else out

def pattern_mdl_compression(hyperedges_df: pd.DataFrame) -> pd.DataFrame:
    if hyperedges_df.empty: return pd.DataFrame()
    df=hyperedges_df.copy(); total=max(df['Count'].sum(),2)
    df['Bits_naive']=df['Size']*np.log2(total); df['Bits_mdl']=-np.log2(df['Support'].clip(lower=1e-12)); df['Compression_Gain']=df['Bits_naive']-df['Bits_mdl']
    return df.sort_values('Compression_Gain', ascending=False)

def shapley_pattern_contributions(class_df: pd.DataFrame, threshold: int=3, top_n:int=10) -> pd.DataFrame:
    prev=class_df.mean(skipna=True).sort_values(ascending=False); feats=prev.head(min(top_n,len(prev))).index.tolist()
    if not feats: return pd.DataFrame()
    p=prev.loc[feats].fillna(0.5).to_dict(); from math import factorial; n=len(feats); rows=[]
    for f in feats:
        phi=0.0; others=[x for x in feats if x!=f]
        for r in range(len(others)+1):
            for S in combinations(others,r):
                probS=sum(p[x] for x in S); vS=1.0 if probS>=threshold else 0.0; vSf=1.0 if (probS+p[f])>=threshold else 0.0
                phi += factorial(r)*factorial(n-r-1)/factorial(n)*(vSf-vS)
        rows.append({'Class':f,'Shapley_MDR_Contribution':phi,'Prevalence':p[f]})
    return pd.DataFrame(rows).sort_values('Shapley_MDR_Contribution', ascending=False)
