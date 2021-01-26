#BCC
#LUCAS SOUTO
#JOSUÉ MELO

from visual import *
from visual.graph import *
import random as rand
from random import random

#VARIAVEIS ESTOCASTICAS
Natoms = rand.randrange(20,300)      # quantidade min e max aleatoria de atomos
L = rand.randrange(1,5)             # medidas do cubo
T = rand.randrange(1000,10000)      # temperatura

print("Digite a simulação que deseja: ")
x = input(
        "Digite 1: Temperatura baixa, nº de atomos baixo e tamanho do recipiente pequeno\n" +
        "Digite 2 Temperatura média, nº de atomos médio e tamanho do recipiente médio\n" +
        "Digite 3 Temperatura alta, nº de atomos alto e tamanho do recipiente grande\n")
            
if x==1:
    T = rand.randrange(10,100)
    L = rand.randrange(1,3)
    Natoms = rand.randrange(10,30)
elif x==2 :
    T = rand.randrange(1000,5000)
    L = rand.randrange(2,5)
    Natoms = rand.randrange(50,100)
elif x==3 :
    T = rand.randrange(10000,50000)
    L = rand.randrange(5,10)
    Natoms = rand.randrange(500,1000)
else:
    print("Entrada para aleatório!")
             
win=700 #tamanho da tela

gray = (0.8,0.8,0.8) # cor das bordas do cubo
Matom = 4E-3/6E23 # massa
Ratom = 0.03 # tamanho exagerado do átomo de hélio
k = 1.4E-23 # constante de Boltzmann (pressão do gás)
dt = 1E-5 #passos tempo

scene = display(title="Gas ideal", width=win, height=win, x=0, y=0, center=(L/2.,L/2.,L/2.))

xaxis = curve(pos=[(0,0,0), (L,0,0)], color=gray)
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=gray)
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=gray)
xaxis2 = curve(pos=[(L,L,L), (0,L,L), (0,0,L), (L,0,L)], color=gray)
yaxis2 = curve(pos=[(L,L,L), (L,0,L), (L,0,0), (L,L,0)], color=gray)
zaxis2 = curve(pos=[(L,L,L), (L,L,0), (0,L,0), (0,L,L)], color=gray)

Atoms = []
poslist = []
plist = []
mlist = []
rlist = []

for i in range(Natoms):
    Lmin = 1.1*Ratom
    Lmax = L-Lmin
    
    x = Lmin+(Lmax-Lmin)*random()
    y = Lmin+(Lmax-Lmin)*random()
    z = Lmin+(Lmax-Lmin)*random()
    r = Ratom
    
    Atoms.append(sphere(pos=(x,y,z), radius=r, color=color.red))
    mass = Matom*r**3/Ratom**3
    
    pavg = sqrt(2.*mass*1.5*k*T) 
    theta = pi*random()
    phi = 2*pi*random()
    
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    
    poslist.append((x,y,z))
    
    plist.append((px,py,pz))
    mlist.append(mass)
    rlist.append(r)

pos = array(poslist)
p = array(plist)
m = array(mlist)
m.shape = (Natoms,1)
radius = array(rlist)

pos = pos+(p/m)*(dt/2.)


print("Variaveis estocasticas:")
print("Numero de atomos: " + str(Natoms))
print("Tamanho do cubo: " + str(L))
print("Temperatura: " + str(T) + "K")


while True:
    rate(50)
    
    # atualiza posições
    pos = pos+(p/m)*dt

    r = pos-pos[:,newaxis] 
    rmag = sqrt(sum(square(r),-1)) 
    hit = less_equal(rmag,radius+radius[:,None])-identity(Natoms)
    hitlist = sort(nonzero(hit.flat)[0]).tolist()

    for ij in hitlist:
        i, j = divmod(ij,Natoms)
        hitlist.remove(j*Natoms+i) 
        ptot = p[i]+p[j]
        mi = m[i,0]
        mj = m[j,0]
        vi = p[i]/mi
        vj = p[j]/mj
        ri = Atoms[i].radius
        rj = Atoms[j].radius
        a = mag(vj-vi)**2
        if a == 0: continue 
        b = 2*dot(pos[i]-pos[j],vj-vi)
        c = mag(pos[i]-pos[j])**2-(ri+rj)**2
        d = b**2-4.*a*c
        if d < 0: continue 
        deltat = (-b+sqrt(d))/(2.*a) 
        pos[i] = pos[i]-(p[i]/mi)*deltat 
        pos[j] = pos[j]-(p[j]/mj)*deltat
        mtot = mi+mj
        pcmi = p[i]-ptot*mi/mtot 
        pcmj = p[j]-ptot*mj/mtot
        rrel = norm(pos[j]-pos[i])
        pcmi = pcmi-2*dot(pcmi,rrel)*rrel 
        pcmj = pcmj-2*dot(pcmj,rrel)*rrel
        p[i] = pcmi+ptot*mi/mtot 
        p[j] = pcmj+ptot*mj/mtot
        pos[i] = pos[i]+(p[i]/mi)*deltat
        pos[j] = pos[j]+(p[j]/mj)*deltat
 
    # Parede
    outside = less_equal(pos,Ratom) 
    p1 = p*outside
    p = p-p1+abs(p1) 
    outside = greater_equal(pos,L-Ratom) 
    p1 = p*outside
    p = p-p1-abs(p1) 

    
    

    for i in range(Natoms):
        Atoms[i].pos = pos[i]


