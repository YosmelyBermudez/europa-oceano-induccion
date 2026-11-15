import numpy as np, pandas as pd
RJ = 71492.0
def leer_sts(path):
    lines = open(path).read().splitlines()
    i = max(k for k,l in enumerate(lines) if l.strip()=='END_OBJECT')
    rows = [l.split() for l in lines[i+1:] if l.strip()]
    a = np.array(rows, float)
    t = pd.to_datetime([f"{int(r[0])}-{int(r[1]):03d}" for r in a], format="%Y-%j") \
        + pd.to_timedelta(a[:,2],'h')+pd.to_timedelta(a[:,3],'m')+pd.to_timedelta(a[:,4]+a[:,5]/1e3,'s')
    df = pd.DataFrame({'t':t,'Bx':a[:,7],'By':a[:,8],'Bz':a[:,9],
                       'x':a[:,11]/RJ,'y':a[:,12]/RJ,'z':a[:,13]/RJ})
    return df
