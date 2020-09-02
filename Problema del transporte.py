# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 22:20:00 2020

@author: Juan
"""

# -*- coding: utf-8 -*-
import numpy as np



#  Creado por Juan Antonio Caro y Luis Cortés

matrizEj = [[1200,8500,1850,2250],[400,800,900,1400],[800,1200,1000,1100]]
listaS = [2500000,1800000,1500000]
listaD = [2000000,1500000,1200000,1200000]
mc = [[2,3,4,5],[5,4,3,1],[1,3,3,2]]
ld1 = [6,11,17,12]

def nwRule(ls,ld):
    
    ls = list(ls)
    ld = list(ld)
    res = [[0]*len(ld) for i in range(len(ls))]
    
    i=0
    j=0
    while (i < len(ls) and j < len(ld)):
        if(ls[i] >= ld[j]):
            res[i][j] = ld[j]
            ls[i] = ls[i] - ld[j]
            j += 1
            if (ls[i] == 0):
                i += 1
        else:
            res[i][j] = ls[i]
            ld[j] = ld[j] - ls[i]
            i += 1

    return res

def muRule(m,ls,ld):

    mc = []
    for i in m:
        mc.append(list(i))
    ls = list(ls)
    ld = list(ld)
    
    def escogeVariable(mc):
        minimo =  float('inf')
        res = 0
        for i in mc:
            if(minimo > min(i)):
                res = (mc.index(i), i.index(min(i)))
                minimo = min(i)
        return res

    res = [[0]*len(ld) for i in range(len(ls))]
    coordenadas = escogeVariable(mc)
    sld = sum(i for i in ld)
    sls = sum(i for i in ls)

    while(sld != 0 and sls != 0):
        i = coordenadas[0]
        j = coordenadas[1]

        if(ls[i] >= ld[j]):
            res[i][j] = ld[j]
            ls[i] = ls[i] - ld[j]
            mc[i][j] = ld[j]
            ld[j] = 0
            for k in range(len(mc)):
                    mc[k][j] = float('inf')
        else:
            res[i][j] = ls[i]
            ld[j] = ld[j] - ls[i]
            mc[i][j] = ls[i]
            ls[i] = 0
            for k in range(len(mc[i])):
                    mc[i][k] = float('inf')

        sld = sum(i for i in ld)
        sls = sum(i for i in ls)
        coordenadas = escogeVariable(mc)

    return res

def vogel(mC,ls,ld):

    def escogeVariable(mc):

        ucd1 = [0]*len(mc)
        for i in range(len(mc)):
                l = sorted(mc[i])
                if(l[1] != float('inf')):
                    ucd1[i] = l[1] - l[0]
                else:
                    if(l[0] != float('inf')):
                        ucd1[i] = 0
                    else:
                        ucd1[i] = -1

        ucd2 = [0]*len(mc[0])
        traspuesta = list(zip(*mc))
        for j in range(len(traspuesta)):
            l = sorted(traspuesta[i])
            if(l[1] != float('inf')):
                ucd2[j] = l[1] - l[0]
            else:
                if(l[0] != float('inf')):
                    ucd2[j] = 0
                else:
                    ucd2[j] = -1

        ucd = ucd1 + ucd2
        maximo = ucd.index(max(ucd))

        if (maximo < len(mc)):
            x = ucd1.index(max(ucd1))
            y = mc[x].index(min(mc[x]))
        else:
            maximo = maximo - len(mc)
            x = traspuesta[maximo].index(min(traspuesta[maximo]))
            y = maximo

        return (x,y)


    ls = list(ls)
    ld = list(ld)
    mc = []
    for x in mC:
        mc.append(list(x))
    res = [[0]*len(ld) for i in range(len(ls))]
    par = escogeVariable(mc)

    sld = sum(i for i in ld)
    sls = sum(i for i in ls)

    while(sld != 0 and sls != 0):
        i = par[0]
        j = par[1]
        if(ls[i] >= ld[j]):
            res[i][j] = ld[j]
            ls[i] = ls[i] - ld[j]
            mc[i][j] = ld[j]
            ld[j] = 0
            for k in range(len(mc)):
                    mc[k][j] = float('inf')
        else:
            res[i][j] = ls[i]
            ld[j] = ld[j] - ls[i]
            mc[i][j] = ls[i]
            ls[i] = 0
            for k in range(len(mc[i])):
                    mc[i][k] = float('inf')
        sld = sum(i for i in ld)
        sls = sum(i for i in ls)
        par = escogeVariable(mc)
        
    return res

def optimizador(mc, ls, ld, sol):

    def calcula_uv(m,s):
        A = []
        b = []
        for i in range(len(m)):
            for j in range (len(m[0])):
                l = [0]*(len(m)+len(m[0]))
                if(s[i][j] != 0):
                    l[j]=1
                    l[i+len(m[0])]=1
                    A.append(l)
                    b.append(m[i][j])
        l = [0]*(len(m)+len(m[0]))
        l[0]=1
        A.append(l)
        b.append(0)
        return list(np.linalg.solve(A,b))

    def resta_uv(m,s):
        res = [[0]*len(m[0]) for i in range(len(m))]
        uv = calcula_uv(m,s)
        for i in range(len(mc)):
            for j in range (len(m[0])):
                res[i][j] = m[i][j] - uv[j] - uv[i+len(m[0])]
        
        return res

    def findminimum(m,s):
        matrix = resta_uv(m,s)
        minimum = min(min(i) for i in matrix)
        for i in range(len(mc)):
            for j in range (len(m[0])):
                if matrix[i][j] == minimum:
                    return (i,j)

    def findpath(m,s):
        
        a,b = findminimum(m,s)
        s[a][b] = 1 # Iniciamos el camino
        m1 = [len(m[0])*[0] for i in range(len(m))]
        m1[a][b] = 1
        camino = [(a,b)]

        def direccion(camino):
            res = 0
            if len(camino) >1:
                p = camino[-1][1] == camino[-2][1]
                if (p == 0):
                    res = 0 # se moverá en vertical
                else:
                    res = 1 # se moverá en horizontal
            return res

        def findpath2(camino):
            d = direccion(camino)
            p = camino[-1]
                
            def incrementa_todo(d,p):
                res = []
                if (d == 0):
                    for x in range(len(m)-1):
                        res.append(((p[0]+x+1) % (len(m)),p[1]))
                else:
                    for x in range(len(m[0])-1):
                        res.append((p[0],(p[1]+x+1) % len(m[0])))
                return res
            
            def camino_valido(d,sp):
                res = False
                if (d == 0):
                    k = sp[0]
                    res = sum([s[k][i]!=0 for i in range(len(s[0]))]) > 1
                else:
                    k = sp[1]
                    res = sum([i[k]!=0 for i in s]) > 1
                return res

            
            opciones = incrementa_todo(d,p)
            opciones = [i for i in opciones if camino_valido(d,i) and s[i[0]][i[1]] != 0]
            res = False
            if not opciones:
                return False
            else:
                if (camino[0] in opciones):
                    return camino
                else:
                    while not res:
                        path = list(camino)
                        c = opciones.pop(0)
                        path.append(c)
                        res = findpath2(path)
                    return res

            
    
        camino = findpath2(camino)
        s[a][b] = 0 #Restaurar a 0 la posición previamente modificada
        return camino
    
    def modificador(lista,minimo):
        alternador = 0
        for a,b in lista:
            sol[a][b] += (1 - alternador) * minimo - (alternador) * minimo
            alternador = (alternador + 1) % 2 
            
    def negativo(m,s):
        res = min([min(x) for x in resta_uv(m,s)])
        return res < 0
    
    while negativo(mc,sol):
        solucionado = findpath(mc,sol)
        minimizar = [solucionado[i] for i in range(len(solucionado)) if i%2 == 1]
        minimo = min([sol[a][b] for a,b in minimizar])
        modificador(solucionado,minimo)
    

    return sol

def costeOptimal(costes,solucion):
    return sum([costes[x][y] * solucion[x][y] for y in range(len(costes[0])) for x in range(len(costes))])


mc = [[2,3,4,5],[5,4,3,1],[1,3,3,2]]
ls = [10,15,21]
ld = [6,11,17,12]

print("Ejecutando test1.py")

print("Matriz de costes:")
for i in range(len(mc)):
    print(mc[i])

print("Lista de origen")
for i in range(len(ls)):
    print([ls[i]])

print("Lista de destino")
print(ld)

print("Resultado de la ejecución de muRule:")
res = muRule(mc,ls,ld)
for i in range(len(res)):
    print(res[i])

print("Resultado de la ejecución de vogel:")
res = vogel(mc,ls,ld)
for i in range(len(res)):
    print(res[i])

print("Resultado de la ejecución de nwRule:")
res = nwRule(ls,ld)
for i in range(len(res)):
    print(res[i])
    
print("Resultado al optimizar nwRule:")
res = optimizador(mc,ls,ld,res)
for i in range(len(res)):
    print(res[i])
    
print("Coste de la solución optimizada:")
res = costeOptimal(mc,res)
print(res)

mc = [[1200,8500,1850,2250],[400,800,900,1400],[800,1200,1000,1100]]
ls = [2500000,1800000,1500000]
ld = [2000000,1500000,1200000,1200000]

print("Ejecutando test2.py")

print("Matriz de costes:")
for i in range(len(mc)):
    print(mc[i])

print("Lista de origen")
for i in range(len(ls)):
    print([ls[i]])

print("Lista de destino")
print(ld)

print("Resultado de la ejecución de muRule:")
res = muRule(mc,ls,ld)
for i in range(len(res)):
    print(res[i])

print("Resultado de la ejecución de vogel")
res = vogel(mc,ls,ld)
for i in range(len(res)):
    print(res[i])

print("Resultado de la ejecución de nwRule:")
res = nwRule(ls,ld)
for i in range(len(res)):
    print(res[i])
    
print("Resultado al optimizar nwRule:")
res = optimizador(mc,ls,ld,res)
for i in range(len(res)):
    print(res[i])
    
print("Coste de la solución optimizada:")
res = costeOptimal(mc,res)
print(res)