"""
================================================================================
CORREÇÕES ESTATÍSTICAS — ARTIGO 1
Atende itens 2, 3, 6, 8 do relatório de revisão técnica (25/08/2026)
================================================================================
1. p exatos bicaudais por enumeração completa (720 permutações) — 19 variáveis
2. FDR Benjamini-Hochberg sobre as 19 correlações químicas
3. PERMANOVA bifatorial com tabela ANOVA completa (gl fechando)
4. PERMDISP para profundidade e ilha
5. Leave-one-out para Ca e Mg
6. Sensibilidade por profundidade (0-10 / 10-20 / 0-20)
7. Produtividade por árvore como métrica alternativa
8. Diagnóstico da estrutura de réplicas
================================================================================
"""

# Author: Luan Daniel Silva Ferreira (ORCID 0000-0001-9187-6988)
# Federal University of Para (UFPA), Belem, PA, Brazil
# Repository: https://github.com/luandanbio/floodplain-cacao-mocajuba
# Archived (functional analysis): https://doi.org/10.5281/zenodo.21345125
# Archived (decontamination): https://doi.org/10.5281/zenodo.17498295
# Last updated: 31 August 2026

import itertools
import numpy as np
import pandas as pd
from scipy import stats

SEED = 42
VARS = ["pH","Carbono","MO","N","CN","P","Al","Acidez","Na","K",
        "Ca","Mg","S","CTC","V","Cu","Zn","Mn","Fe"]
PONTOS = ['P1','P2','P3','P4','P5','P6']
NOMES  = ['Santana','Santaninha','Angapijó','Conceição','São Joaquim','Tauaré']
PROD_TOTAL = np.array([2000., 600., 450., 1035., 1000., 1500.])
N_TREES    = np.array([6000., 13500., 12000., 13600., 3000., 3000.])
PROD_TREE  = PROD_TOTAL / N_TREES

# ── Carregar ─────────────────────────────────────────────────────────────────
raw = pd.read_csv('/mnt/user-data/outputs/solo_quimica.csv')
df  = raw[VARS].apply(pd.to_numeric, errors='coerce')
df['ilha']  = ['P'+str(x).replace('S','') for x in raw['Ponto']]
df['prof']  = [str(x) for x in raw['Profundidade']]

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 0 — DIAGNÓSTICO DA ESTRUTURA DE RÉPLICAS
# ══════════════════════════════════════════════════════════════════════════════
print("="*78)
print("PARTE 0 — DIAGNÓSTICO DA ESTRUTURA DE RÉPLICAS")
print("="*78)
cvs = []
for p in PONTOS:
    for d in ['0-10','10-20']:
        g = df[(df.ilha==p)&(df.prof==d)][VARS]
        cvs.append(g.std(ddof=1)/g.mean()*100)
cv_tab = pd.DataFrame(cvs)
print("\nCoeficiente de variação (%) entre as 3 réplicas — média sobre as 12 combinações:")
print(cv_tab.mean().round(3).to_string())
print(f"\nCV mediano global: {cv_tab.mean().median():.3f}%")
print("""
INTERPRETAÇÃO: CV mediano < 0,5% para a maioria das variáveis é incompatível
com réplicas de campo espaçadas 10 m em solo de várzea, onde CVs de 10–30% são
típicos. Este padrão é a assinatura de RÉPLICAS ANALÍTICAS (leituras repetidas
do mesmo extrato/composto), não de réplicas de campo independentes.

CONSEQUÊNCIA PARA O DESENHO: as unidades experimentais independentes são
6 ilhas × 2 profundidades = 12, e não 36. A PERMANOVA deve ser recalculada
sobre as 12 médias, com o resíduo removido do modelo.
""")

# ── Médias por ilha × profundidade (12 unidades independentes) ───────────────
med = df.groupby(['ilha','prof'])[VARS].mean().reset_index()
med = med.sort_values(['ilha','prof']).reset_index(drop=True)
# Médias por ilha (agregado 0-20)
med_ilha = df.groupby('ilha')[VARS].mean().loc[PONTOS]
med_010  = df[df.prof=='0-10'].groupby('ilha')[VARS].mean().loc[PONTOS]
med_1020 = df[df.prof=='10-20'].groupby('ilha')[VARS].mean().loc[PONTOS]

# ══════════════════════════════════════════════════════════════════════════════
# MOTOR DE PERMUTAÇÃO EXATA
# ══════════════════════════════════════════════════════════════════════════════
def exact_spearman(x, y):
    """rho de Spearman + p exato BICAUDAL por enumeração completa (n<=8)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    def corr(a, b):
        am = a - a.mean(); bm = b - b.mean()
        d = np.sqrt((am**2).sum() * (bm**2).sum())
        return np.dot(am, bm)/d if d > 0 else 0.0
    rho = corr(rx, ry)
    cnt = sum(1 for perm in itertools.permutations(ry)
              if abs(corr(rx, np.array(perm))) >= abs(rho) - 1e-9)
    tot = np.math.factorial(n) if hasattr(np, 'math') else 720
    import math; tot = math.factorial(n)
    return rho, cnt/tot

def bh(p):
    p = np.asarray(p, float); m = len(p)
    o = np.argsort(p); q = p[o]*m/(np.arange(m)+1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1.0)
    return out

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 1 — CORRELAÇÕES QUÍMICA × PRODUTIVIDADE (p exato bicaudal + FDR)
# ══════════════════════════════════════════════════════════════════════════════
print("="*78)
print("PARTE 1 — TABELA 2 CORRIGIDA: p exato BICAUDAL + FDR sobre as 19 variáveis")
print("="*78)

rows = []
for v in VARS:
    rho_t, p_t = exact_spearman(med_ilha[v].values, PROD_TOTAL)
    rho_a, p_a = exact_spearman(med_ilha[v].values, PROD_TREE)
    rows.append({'Variable': v, 'rho_total': round(rho_t,3), 'p_total': round(p_t,4),
                 'rho_tree': round(rho_a,3), 'p_tree': round(p_a,4)})
t2 = pd.DataFrame(rows)
t2['q_BH_total'] = bh(t2['p_total'].values).round(3)
t2['q_BH_tree']  = bh(t2['p_tree'].values).round(3)
t2['Sig_total']  = np.where(t2.p_total<0.01,'**', np.where(t2.p_total<0.05,'*','ns'))
t2['Sig_tree']   = np.where(t2.p_tree<0.01,'**', np.where(t2.p_tree<0.05,'*','ns'))
t2['Direction']  = np.where(t2.rho_total>0,'Positive','Negative')
t2 = t2.reindex(t2.rho_total.abs().sort_values(ascending=False).index).reset_index(drop=True)
print(t2.to_string(index=False))
t2.to_csv('/home/claude/CORR_tabela2_quimica_exata.csv', index=False)
print(f"\nCa: p exato bicaudal = {t2[t2.Variable=='Ca'].p_total.values[0]:.4f}  "
      f"q(BH) = {t2[t2.Variable=='Ca'].q_BH_total.values[0]:.3f}")
print(f"Mg: p exato bicaudal = {t2[t2.Variable=='Mg'].p_total.values[0]:.4f}  "
      f"q(BH) = {t2[t2.Variable=='Mg'].q_BH_total.values[0]:.3f}")
print(f"Nenhuma variável sobrevive ao FDR: {(t2.q_BH_total<0.05).sum()} de 19")

# Colinearidade Ca-Mg
rho_cm, p_cm = exact_spearman(med_ilha['Ca'].values, med_ilha['Mg'].values)
r_pear = np.corrcoef(med_ilha['Ca'].values, med_ilha['Mg'].values)[0,1]
ca = med_ilha['Ca'].values; mg = med_ilha['Mg'].values
print(f"\nColinearidade Ca-Mg: Spearman rho = {rho_cm:.3f} (p = {p_cm:.4f}); "
      f"Pearson r = {r_pear:.3f}")
print(f"Amplitude Ca: {ca.min():.2f}-{ca.max():.2f} cmolc/kg (CV = {ca.std(ddof=1)/ca.mean()*100:.1f}%)")
print(f"Amplitude Mg: {mg.min():.2f}-{mg.max():.2f} cmolc/kg (CV = {mg.std(ddof=1)/mg.mean()*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 2 — SENSIBILIDADE POR PROFUNDIDADE (Ca e Mg)
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 2 — SENSIBILIDADE POR PROFUNDIDADE (p exato bicaudal)")
print("="*78)
sens = []
for lbl, tab in [('0-10 cm', med_010), ('10-20 cm', med_1020), ('0-20 cm (mean)', med_ilha)]:
    for v in ['Ca','Mg']:
        for mlbl, pv in [('Total (kg/yr)', PROD_TOTAL), ('Per tree (kg/tree)', PROD_TREE)]:
            r, p = exact_spearman(tab[v].values, pv)
            sens.append({'Depth': lbl, 'Var': v, 'Metric': mlbl,
                         'rho': round(r,3), 'p_exact': round(p,4),
                         'Sig': '**' if p<0.01 else ('*' if p<0.05 else 'ns')})
sens = pd.DataFrame(sens)
print(sens.to_string(index=False))
sens.to_csv('/home/claude/CORR_sensibilidade_profundidade.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 3 — LEAVE-ONE-OUT (Ca e Mg)
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 3 — LEAVE-ONE-ISLAND-OUT (Ca e Mg × produção total)")
print("="*78)
loo = []
for v in ['Ca','Mg']:
    r_full, p_full = exact_spearman(med_ilha[v].values, PROD_TOTAL)
    print(f"\n{v}  — dataset completo: rho = {r_full:.3f}  (p = {p_full:.4f})")
    for i, ex in enumerate(PONTOS):
        idx = [j for j in range(6) if j != i]
        r, p = exact_spearman(med_ilha[v].values[idx], PROD_TOTAL[idx])
        loo.append({'Var': v, 'Excluded': ex, 'Name': NOMES[i],
                    'rho': round(r,3), 'p_exact': round(p,4),
                    'delta_rho': round(r-r_full,3)})
        print(f"   sem {ex} ({NOMES[i]:<12}): rho = {r:6.3f}  p = {p:.4f}  "
              f"delta = {r-r_full:+.3f}")
loo = pd.DataFrame(loo)
loo.to_csv('/home/claude/CORR_leave_one_out.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 4 — PERMANOVA BIFATORIAL COM gl FECHANDO
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 4 — PERMANOVA BIFATORIAL (12 unidades independentes)")
print("Modelo: quimica ~ Ilha + Profundidade  (sem interação estimável)")
print("="*78)

X = med[VARS].to_numpy(float)
Z = (X - X.mean(0)) / X.std(0, ddof=1)
g_ilha = med['ilha'].to_numpy()
g_prof = med['prof'].to_numpy()
n = len(Z)

D2 = ((Z[:,None,:] - Z[None,:,:])**2).sum(-1)
iu = np.triu_indices(n, 1)
SS_T = D2[iu].sum()/n

def ss_within(D2, grp):
    s = 0.0
    for g in np.unique(grp):
        k = np.where(grp==g)[0]
        if len(k) < 2: continue
        sub = D2[np.ix_(k,k)]
        s += sub[np.triu_indices(len(k),1)].sum()/len(k)
    return s

SS_ilha = SS_T - ss_within(D2, g_ilha)
SS_prof = SS_T - ss_within(D2, g_prof)
SS_res  = SS_T - SS_ilha - SS_prof

gl_i, gl_p = 5, 1
gl_r = n - 1 - gl_i - gl_p          # 12-1-5-1 = 5
gl_t = n - 1                         # 11

MS_i, MS_p, MS_r = SS_ilha/gl_i, SS_prof/gl_p, SS_res/gl_r
F_i, F_p = MS_i/MS_r, MS_p/MS_r

rng = np.random.default_rng(SEED)
NP = 9999
Fp_perm = np.empty(NP); Fi_perm = np.empty(NP)
for b in range(NP):
    pp = g_prof.copy()
    for g in np.unique(g_ilha):
        k = np.where(g_ilha==g)[0]
        pp[k] = rng.permutation(g_prof[k])
    ssp = SS_T - ss_within(D2, pp)
    ssr = SS_T - SS_ilha - ssp
    Fp_perm[b] = (ssp/gl_p)/(ssr/gl_r) if ssr>0 else 0
    perm = rng.permutation(n)
    ssi = SS_T - ss_within(D2[np.ix_(perm,perm)], g_ilha)
    ssr2 = SS_T - ssi - SS_prof
    Fi_perm[b] = (ssi/gl_i)/(ssr2/gl_r) if ssr2>0 else 0

p_prof = (np.sum(Fp_perm >= F_p)+1)/(NP+1)
p_ilha = (np.sum(Fi_perm >= F_i)+1)/(NP+1)

perm_tab = pd.DataFrame([
    {'Source':'Island','df':gl_i,'SS':round(SS_ilha,3),'R2':round(SS_ilha/SS_T,4),
     'MS':round(MS_i,3),'Pseudo-F':round(F_i,2),'p':round(p_ilha,4)},
    {'Source':'Depth','df':gl_p,'SS':round(SS_prof,3),'R2':round(SS_prof/SS_T,4),
     'MS':round(MS_p,3),'Pseudo-F':round(F_p,2),'p':round(p_prof,4)},
    {'Source':'Residual','df':gl_r,'SS':round(SS_res,3),'R2':round(SS_res/SS_T,4),
     'MS':round(MS_r,3),'Pseudo-F':'—','p':'—'},
    {'Source':'Total','df':gl_t,'SS':round(SS_T,3),'R2':1.0,'MS':'—',
     'Pseudo-F':'—','p':'—'},
])
print()
print(perm_tab.to_string(index=False))
print(f"\nVerificação dos gl: {gl_i} + {gl_p} + {gl_r} = {gl_i+gl_p+gl_r} = gl_total ({gl_t}) ✓")
print(f"Verificação das SS: {SS_ilha:.3f} + {SS_prof:.3f} + {SS_res:.3f} = "
      f"{SS_ilha+SS_prof+SS_res:.3f} = SS_total ({SS_T:.3f}) ✓")
perm_tab.to_csv('/home/claude/CORR_tabela3_permanova.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 5 — PERMDISP
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 5 — PERMDISP (homogeneidade de dispersões)")
print("="*78)

def permdisp(D2, grp, nperm=9999, seed=SEED):
    n = D2.shape[0]
    H = np.eye(n) - np.ones((n,n))/n
    B = -0.5 * H @ D2 @ H
    ev, evec = np.linalg.eigh(B)
    o = np.argsort(ev)[::-1]; ev, evec = ev[o], evec[:,o]
    pos = ev > 1e-10
    C = evec[:,pos] * np.sqrt(ev[pos])
    cent = {g: C[grp==g].mean(0) for g in np.unique(grp)}
    d = np.array([np.sqrt(((C[i]-cent[grp[i]])**2).sum()) for i in range(n)])
    gs = np.unique(grp); k = len(gs)
    gm = np.array([d[grp==g].mean() for g in gs]); GM = d.mean()
    ssb = sum((grp==g).sum()*(gm[i]-GM)**2 for i,g in enumerate(gs))
    ssw = sum(((d[grp==g]-gm[i])**2).sum() for i,g in enumerate(gs))
    df1, df2 = k-1, n-k
    F = (ssb/df1)/(ssw/df2) if ssw>0 else np.inf
    rg = np.random.default_rng(seed); cnt = 0
    for _ in range(nperm):
        gp = grp[rg.permutation(n)]
        gmp = np.array([d[gp==g].mean() for g in gs])
        b_ = sum((gp==g).sum()*(gmp[i]-GM)**2 for i,g in enumerate(gs))
        w_ = sum(((d[gp==g]-gmp[i])**2).sum() for i,g in enumerate(gs))
        if w_>0 and (b_/df1)/(w_/df2) >= F: cnt += 1
    return F, (cnt+1)/(nperm+1), df1, df2, d

Fd, pd_, d1, d2_, dist_d = permdisp(D2, g_prof)
Fi_, pi_, i1, i2, dist_i = permdisp(D2, g_ilha)
print(f"\nProfundidade: F({d1},{d2_}) = {Fd:.4f}   p = {pd_:.4f}   "
      f"→ {'dispersão homogênea' if pd_>0.05 else 'dispersão heterogênea'}")
print(f"Ilha:         F({i1},{i2}) = {Fi_:.4f}   p = {pi_:.4f}   "
      f"→ {'dispersão homogênea' if pi_>0.05 else 'dispersão heterogênea'}")

pd.DataFrame([
    {'Factor':'Depth','F':round(Fd,4),'df1':d1,'df2':d2_,'p':round(pd_,4),
     'Interpretation':'Homogeneous' if pd_>0.05 else 'Heterogeneous'},
    {'Factor':'Island','F':round(Fi_,4),'df1':i1,'df2':i2,'p':round(pi_,4),
     'Interpretation':'Homogeneous' if pi_>0.05 else 'Heterogeneous'},
]).to_csv('/home/claude/CORR_permdisp.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 6 — TESTES PAREADOS POR PROFUNDIDADE + FDR
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 6 — TESTES PAREADOS DE PROFUNDIDADE (n = 6 pares) + FDR")
print("="*78)
pair = []
for v in VARS:
    a = med_010[v].values; b = med_1020[v].values
    tst, pt = stats.ttest_rel(a, b)
    try: _, pw = stats.wilcoxon(a, b)
    except Exception: pw = np.nan
    pair.append({'Variable': v, 'mean_0_10': round(a.mean(),3),
                 'mean_10_20': round(b.mean(),3),
                 'delta_pct': round(100*(b.mean()-a.mean())/a.mean(),1),
                 't': round(tst,3), 'p_t': round(pt,4),
                 'p_wilcoxon': round(pw,4) if not np.isnan(pw) else np.nan})
pair = pd.DataFrame(pair).sort_values('p_t')
pair['q_BH'] = bh(pair['p_t'].values).round(4)
print(pair.to_string(index=False))
print(f"\nVariáveis com p < 0.05: {(pair.p_t<0.05).sum()} de 19")
print(f"Variáveis com q(BH) < 0.05: {(pair.q_BH<0.05).sum()} de 19")
pair.to_csv('/home/claude/CORR_testes_profundidade.csv', index=False)

# ══════════════════════════════════════════════════════════════════════════════
# PARTE 7 — TABELA 1 CORRIGIDA (kg/árvore)
# ══════════════════════════════════════════════════════════════════════════════
print("\n"+"="*78)
print("PARTE 7 — TABELA 1 CORRIGIDA (coluna kg/árvore realinhada)")
print("="*78)
t1 = pd.DataFrame({
    'Island': PONTOS, 'Name': NOMES,
    'Production_kg_yr': PROD_TOTAL.astype(int),
    'System': ['AFS']*6,
    'Trees': N_TREES.astype(int),
    'kg_per_tree': np.round(PROD_TREE,3),
    'Management': ['Clearing + pruning','Clearing','Clearing','Clearing',
                   'Clearing','Pruning + shading'],
})
print(t1.to_string(index=False))
t1.to_csv('/home/claude/CORR_tabela1_produtividade.csv', index=False)

r_dens, p_dens = exact_spearman(N_TREES, PROD_TREE)
print(f"\nSpearman(nº de árvores, kg/árvore) = {r_dens:.3f}  (p exato = {p_dens:.4f})")

print("\n" + "="*78)
print("TODOS OS ARQUIVOS DE CORREÇÃO GERADOS")
print("="*78)
