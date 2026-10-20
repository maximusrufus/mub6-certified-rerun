import sys
sys.argv=['x']
from verify_excluded import *
import math
lines=[l.split() for l in open('../stageA/tests/sample_seeds.txt')]
rows=[l.split() for l in open('prof/ex.txt') if l[0]!='#']
for idx in (0,5):
    q=[((float(b)+0.5)/112.0) for b in lines[idx][1:10]]
    S=Seed(q)
    ex=[r for r in rows if int(r[0])==idx][0]
    Nu=int(ex[2]); Rs=float.fromhex(ex[10]); z2=float.fromhex(ex[8])
    h,z2r,Er=S.level_bounds(Nu,Rs,z2)
    print(idx,Nu,'claim z2',z2,'mine',z2r.upper(),'Rs',Rs,float(S.Rs_true.upper().mid().str(20,radius=False)))
    print(' m',S.m.lower(),S.mB.lower(),S.mC.lower())
def F(x): return float(x.mid().str(20,radius=False))
A0=F(S.A0); mC=F(S.mC); m=min(F(S.mB),F(S.mC)); Rs=float.fromhex(ex[10]); h=math.pi/Nu
Ast=min(2,A0+Rs); qq=Rs**2*(1+(Ast**2+6)/mC**2); p=qq/(math.sqrt(6)+math.sqrt(6-qq)); off=Rs*(math.sqrt(6)+Ast)/m
print('A0',A0,'Ast',Ast,'qq',qq,'p',p,'z2 float',off**2+(6*h+p)**2,'claim',z2)
Rb=arb(Rs); t=S.A0+Rb
Astar = arb.union(min_ball(t.lower(), arb(2)), min_ball(t.upper(), arb(2)))
print('Astar',Astar, 'm',S.m)
qqb=Rb**2*(1+(Astar**2+6)/S.mC**2); print('qq',qqb)
off=Rb*(S.s6+Astar)/S.m; print('off',off, 'float',off)
