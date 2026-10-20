import sys
sys.argv=['x']
from verify_excluded import *
q=[((float(b)+0.5)/112.0) for b in open('../stageA/tests/sample_seeds.txt').readline().split()[1:10]]
S=Seed(q)
ex=[l.split() for l in open('prof/ex.txt') if l[0]!='#'][0]
Nu=int(ex[2]); Rs=float.fromhex(ex[10]); z2=float.fromhex(ex[8])
print(ex[:3], Rs, z2)
h,z2r,Er=S.level_bounds(Nu,Rs,z2)
print(z2r, Er, float.fromhex(ex[9]), float.fromhex(ex[7]))
print(S.m, S.mC, S.A0, S.Rs_true)
print(float.fromhex(ex[9]), Er.upper())
import math
Rs=float.fromhex(ex[10]); z2c=float.fromhex(ex[8]); h=math.pi/Nu
A0=float(S.A0.mid().str(20,radius=False)); m=float(S.m.mid().str(20,radius=False)); mC=float(S.mC.mid().str(20,radius=False)); mB=float(S.mB.mid().str(20,radius=False)); Cop=float(S.Cop.mid().str(20,radius=False))
Ast=min(2,A0+Rs)
qq=Rs**2*(1+(Ast**2+6)/mC**2); p=qq/(math.sqrt(6)+math.sqrt(6-qq))
L=max(math.sqrt(6),Cop*math.sqrt(1+(A0/mB)**2))
eoff=Rs**2*(L+A0+1)/(2*m)+math.sqrt(z2c)*Rs/m
e3=math.sqrt(252)*.5*h*h+p
print('float Er',math.sqrt(eoff**2+e3**2),'claim',float.fromhex(ex[9]))
