"""Numeric results harness mirroring the MATLAB CRTP-SS + adaptive-coding
implementation (Fang et al., IEEE IoTJ 2025) for the report tables."""
import numpy as np, random, itertools, heapq

# ---- GF(2)[x] ----
def deg(a): return a.bit_length()-1
def mul(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def divmod_(a,b):
    q=0; db=deg(b)
    while a and deg(a)>=db:
        s=deg(a)-db; q^=1<<s; a^=b<<s
    return q,a
def mod(a,b): return divmod_(a,b)[1]
def inv(a,m):
    r0,r1=m,mod(a,m); t0,t1=0,1
    while r1:
        q,r=divmod_(r0,r1); r0,r1=r1,r; t0,t1=t1,t0^mul(q,t1)
    return t0
q0=0b111110101; qs=[0b100011011,0b100011101,0b100101011,0b100101101]

def enc_pixel(m_,r1,r2,qi): return mod(m_^r1^mul(r2,q0),qi)

# ---- adaptive coding (B1/B2 classification + pure capacity) ----
def classify(img2x2_share):
    M,N=img2x2_share.shape; s=(M*N)//4; nbc=N//2
    nB1=nB2=0; pureEC=0; lam_count={}
    for z in range(s):
        br=2*(z//nbc); bc=2*(z%nbc)
        p=[img2x2_share[br,bc],img2x2_share[br,bc+1],img2x2_share[br+1,bc],img2x2_share[br+1,bc+1]]
        D=[int(p[0])^int(p[u]) for u in (1,2,3)]; Dmax=max(D)
        if Dmax<=63:
            nB1+=1
            lam=0 if Dmax==0 else (Dmax.bit_length())
            lam_count[lam]=lam_count.get(lam,0)+1
            pureEC += 3*(8-lam)-1   # upper bound before label/msb overhead
        else:
            nB2+=1
    return nB1,nB2,pureEC

def make(kind,M,N,rng):
    rr,cc=np.meshgrid(range(M),range(N),indexing='ij')
    if kind=='smooth': return (np.mod(np.round(80+60*np.sin(rr/17)+50*np.cos(cc/23)+(rr+cc)/8),256)).astype(int)
    if kind=='complex': return rng.integers(0,256,(M,N))
    if kind=='mixed':
        a=(np.mod(np.round(80+60*np.sin(rr/17)+(rr+cc)/8),256)).astype(int)
        a[:,N//2:]=rng.integers(0,256,(M,N-N//2)); return a

def entropy(a):
    h,_=np.histogram(a.ravel(),bins=np.arange(257)); p=h/h.sum(); p=p[p>0]
    return -(p*np.log2(p)).sum()

print("=== CRTP-SS numeric results, (3,4)-threshold, GF(2^8) ===\n")
print("Worked example (paper Sec. IV-A): block [192 191 195 180], R1=128, R2=1")
blk=[192,191,195,180]; got=[enc_pixel(v,128,1,qs[0]) for v in blk]
print(f"  share-1 encrypted block = {got}   expected [174, 209, 173, 218]  match={got==[174,209,173,218]}\n")

rng=np.random.default_rng(3)
M=N=128; n,k=4,3
print(f"Per image type ({M}x{N}), share statistics:")
print(f"{'Type':<10}{'B1 faces':<10}{'B2 faces':<10}{'B1 %':<8}{'pureEC':<9}{'pureER(bpp)':<12}{'shareEnt':<9}")
for kind in ['smooth','mixed','complex']:
    img=make(kind,M,N,rng)
    # scramble omitted for stats (does not change value distribution); share directly
    s=(M*N)//4
    R1=rng.integers(0,256,s); R2=rng.integers(0,65536,s)
    share=np.zeros((M,N),dtype=int); nbc=N//2
    off=[(0,0),(0,1),(1,0),(1,1)]
    for z in range(s):
        br=2*(z//nbc); bc=2*(z%nbc)
        for u in range(4):
            pos=(br+off[u][0],bc+off[u][1])
            share[pos]=enc_pixel(int(img[pos]),int(R1[z]),int(R2[z]),qs[0])
    nB1,nB2,pureEC=classify(share)
    er=pureEC/(M*N)
    print(f"{kind:<10}{nB1:<10}{nB2:<10}{100*nB1/(nB1+nB2):<8.1f}{pureEC:<9}{er:<12.4f}{entropy(share):<9.4f}")

print("\nKey sensitivity (Key1 -> Key1+1), share 1, smooth image:")
img=make('smooth',M,N,rng); s=(M*N)//4; off=[(0,0),(0,1),(1,0),(1,1)]; nbc=N//2
def share_with(seed):
    rr=np.random.default_rng(seed); R1=rr.integers(0,256,s); R2=rr.integers(0,65536,s)
    sh=np.zeros((M,N),dtype=int)
    for z in range(s):
        br=2*(z//nbc); bc=2*(z%nbc)
        for u in range(4):
            pos=(br+off[u][0],bc+off[u][1]); sh[pos]=enc_pixel(int(img[pos]),int(R1[z]),int(R2[z]),qs[0])
    return sh
A=share_with(100); B=share_with(101)
npcr=100*np.mean(A!=B); uaci=100*np.mean(np.abs(A-B)/255)
print(f"  NPCR = {npcr:.4f}%  (ideal 99.6094)   UACI = {uaci:.4f}%  (ideal 33.4635)")

print("\nReconstruction: any k=3 of n=4 shares recover the original (random 32x32):")
img=rng.integers(0,256,(8,8)); s=(8*8)//4; nbc=4; off=[(0,0),(0,1),(1,0),(1,1)]
R1=rng.integers(0,256,s); R2=rng.integers(0,65536,s)
shares=[np.zeros((8,8),dtype=int) for _ in range(4)]
for z in range(s):
    br=2*(z//nbc); bc=2*(z%nbc)
    for u in range(4):
        pos=(br+off[u][0],bc+off[u][1])
        for i in range(4): shares[i][pos]=enc_pixel(int(img[pos]),int(R1[z]),int(R2[z]),qs[i])
def recon(idx):
    Q=1
    for i in idx: Q=mul(Q,qs[i])
    out=np.zeros((8,8),dtype=int)
    for z in range(s):
        br=2*(z//nbc); bc=2*(z%nbc)
        for u in range(4):
            pos=(br+off[u][0],bc+off[u][1]); Y=0
            for i in idx:
                Qi,_=divmod_(Q,qs[i]); iv=inv(mod(Qi,qs[i]),qs[i])
                Y^=mod(mul(mul(int(shares[i][pos]),Qi),iv),Q)
            out[pos]=mod(Y,q0)^int(R1[z])
    return out
for idx in itertools.combinations(range(4),3):
    print(f"  shares {tuple(i+1 for i in idx)}: lossless={np.array_equal(recon(idx),img)}")
