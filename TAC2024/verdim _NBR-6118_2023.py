# Programa para a verificação e o dimensionamento de seções poligonais de concreto armado
"""
Programa para a verificação e o dimensionamento de seções poligonais de concreto armado
Integração das tensões no concreto através do Teorema de Green
Baseado em:
CAMPOS FILHO, A. Dimensionamento e verificação de seções poligonais de concreto armado submetidas à flexo compressão oblíqua. 
Porto Alegre, PPGEC-UFRGS, 2014. 
"""
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define o sistema de equações para o dimensionamento
def dimen(vars, nv, xp, yp, nrc, fcd, il, nb, e , lam, xb, yb, fyd, perc, na, max, may,epsc2, epscu, a1, a2):

    x, alph, astotal = vars
    epss = 0.00
    epsi = 0.00
    
    # calcula os esforços resistentes para o terno {x, alpha, As} inicial
    nr, mrx, mry, epss, epsi = esfor(e, fyd, fcd, nv, xp, yp, nb, xb, yb, perc, x, alph, astotal, b, c, epss, epsi, nrc, il, epsc2, epscu, a1, a2)
    eq1 = max - lam * mrx
    eq2 = may - lam * mry
    eq3 = na - lam * nr
    return [eq1, eq2, eq3]


# Define o sistema de equações para a verificação
def verif(vars, nv, xp, yp, nrc, fcd, il, nb, e , astotal, xb, yb, fyd, perc, na, max, may,epsc2, epscu, a1, a2):
    
    x, alph, lam = vars
    epss = 0.00
    epsi = 0.00
    
    # calcula os esforços resistentes para o terno {x, alpha, As} inicial
    nr, mrx, mry, epss, epsi = esfor(e, fyd, fcd, nv, xp, yp, nb, xb, yb, perc, x, alph, astotal, b, c, epss, epsi, nrc, il, epsc2, epscu, a1, a2)
    eq1 = max - lam * mrx
    eq2 = may - lam * mry
    eq3 = na - lam * nr
    return [eq1, eq2, eq3]


def initial(jox, joy, joxy, lx, ly, na, max, may, fcd, fyd, astotal):
    """
    Inicializa os valores para a solução não-linear do sistema

    """
    pi = math.pi
    pi2 = pi/2.00
    alfa0 = 0.00
    # estimativa inicial do valor da inclinação da linha neutra alpha
    if jox == joy and abs(joxy) > 1e-5:
        alfa0 = pi2
    else:
        if jox != joy:
            alfa0 = math.atan(-2 * joxy / (jox - joy)) / 2.0

    # cálculo dos momentos principais de inércia pjx e pjy
    ca0 = math.cos(alfa0)
    sa0 = math.sin(alfa0)
    ca0 = ca0 * ca0
    sa0 = sa0 * sa0
    ss = joxy * math.sin(2 * alfa0)
    pjx = jox * ca0 + joy * sa0 - ss
    pjy = joy * ca0 + jox * sa0 + ss
    alph = 0.0

    if max == 0:
        alph = -math.copysign(pi2, may)
    else:
        alph = math.atan(may * pjx / (max * pjy))
        if max > 0:
            alph = alph + pi

    alfa0 = alph + alfa0
    
    # estimativa inicial da profundidade da linha neutra x
    
    #uu = 0.50
    #uu = (pi + uu) ** 5
    #uu = uu - int(uu) # gera um número randômico entre 0 e 1
    uu = np.random.uniform(0.00,1.00)
    x = (lx + ly) * uu
    #if np.abs(alfa0)<1e-6:
    #    x = 0.50 * ly
    

    # estimativa inicial da armadura total astotal
    if op == 1:
        as1 = abs(max) / (0.4 * ly * fyd[0]) # dimensionamento simplificado à flexão em torno de x
        as2 = abs(may) / (0.4 * lx * fyd[0]) # dimensionamento simplificado à flexão em torno de y

        if na > 0:
            as3 = na / fyd[0] # dimensionamento para a tração axial
        else:
            as3 = np.maximum(0.0, (np.abs(na) - 0.85 * fcd[0] * area) / fyd[0]) # dimensionamento para a compressão centrada

        astotal = as1 + as2 + as3 # área de aço total

    # estimativa inicial do ângulo de inclinação da linha neutra alpha
    alph = alfa0

    return x, alph, astotal


def esfor(e, fyd, fc, nv, xp, yp, nb, xb, yb, perc, x, alfa, astotal,b, c, epss, epsi, nrc, il, epsc2, epscu, a1, a2):
    """
    Calcula os esforços resistentes {nr, mrx, mry}

    """
    nrz = 0.00
    mrx = 0.00
    mry = 0.00
    nrzti = 0.0
    nrzt = 0.0
    mrks = 0.0
    mret = 0.0
    ks1i = 0.0
    ks2i = 0.0
    ks1ii = 0.0
    ks2ii = 0.0
    
    ksp = np.zeros(nv)
    etp = np.zeros(nv)
    ksb = np.zeros(nb)
    etb = np.zeros(nb)
    epsp = np.zeros(nv)

    ca = math.cos(alfa)
    sa = math.sin(alfa)
    tetsc = 0.0
    tetic = 0.0
    for j1 in range(nv):
        ksp[j1] = xp[j1] * ca + yp[j1] * sa
        etp[j1] = -xp[j1] * sa + yp[j1] * ca
        if etp[j1] > tetsc:
            tetsc = etp[j1]
            tkssc = ksp[j1]
        if etp[j1] < tetic:
            tetic = etp[j1]
            tksic = ksp[j1]

    tetsa = 0.0
    tetia = 0.0
    for j1 in range(nb):
        ksb[j1] = xb[j1] * ca + yb[j1] * sa
        etb[j1] = -xb[j1] * sa + yb[j1] * ca
        if etb[j1] > tetsa:
            tetsa = etb[j1]
            tkssa = ksb[j1]
        if etb[j1] < tetia:
            tetia = etb[j1]
            tksia = ksb[j1]

    # verificação dos domínios de deformação
    h = tetsc - tetic
    dd = tetsc - tetia
    x23 = epscu[0] / (0.01 + epscu[0]) * dd # limite entre os domínios 2 e 3

    if x < x23:
        # domínios 1 e 2
        b = -0.01 / (dd - x)
        c = b * (x - tetsc)
        blx = -0.01 / (dd - x) ** 2
        clx = (dd - tetsc) * blx
        epsi = 0.01
        epss = -0.01 * x / (dd - x)
    else:
        if x < h:
            # domínios 3, 4 e 4a
            b = -epscu[0] / x
            c = -epscu[0] - b * tetsc
            blx = epscu[0] / x ** 2
            clx = -tetsc * blx
            epss = -epscu[0]
            epsi = np.maximum(epscu[0] * (dd - x) / x, 0.00)
        else:
            if x > 1e150:
                # reta b
                b = 0.0
                c = -epsc2[0]
                blx = 1e-100
                clx = 1e-100
                epss = -epsc2[0]
                epsi = -epsc2[0]
                return nrz, mrx, mry
            # domínio 5
            b = -epsc2[0] / (x - (epscu[0] - epsc2[0]) / epscu[0] * h)
            c = b * (x - tetsc)
            blx = epsc2[0] / (x - (epscu[0] - epsc2[0]) / epscu[0] * h) ** 2
            clx = (tetsc - (epscu[0] - epsc2[0]) / epscu[0] * h) * blx
            epss = -epsc2[0] * x / (x - (epscu[0] - epsc2[0]) / epscu[0] * h)
            epsi = -epsc2[0] * (x - h) / (x - (epscu[0] - epsc2[0]) / epscu[0] * h)

    # calcula as deformações em cada vértice da poligonal de concreto
    epsp = b * etp + c

    nrzt = 0.0
    mrks = 0.0
    mret = 0.0
    
    
    # calculo dos esforços resistentes e da matriz de rigidez para as barras de aço
    for j1 in range(nb):
        epsb = b * etb[j1] + c
        fyd_j1 = fyd[j1]
        sig, et = aco(e, epsb, fyd_j1)
        dum1 = perc[j1] * astotal * et * (blx * etb[j1] + clx)
        dum2 = perc[j1] * sig
        nrzti = astotal * dum2
        # esforços resistentes da armadura
        nrzt += nrzti
        mrks += nrzti * etb[j1]
        mret -= nrzti * ksb[j1]
            
    if abs(epss - epsi) <= 1e-10:
        
        # solução para o caso de compressão ou tração centrada
        if epss >= 0:
        # tração centrada: não há contribuição do concreto
            nrz = nrzt
            mrx = mrks * ca - mret * sa
            mry = mrks * sa + mret * ca
            return nrz, mrx, mry, epss, epsi
        # compressão centrada 
        for j1 in range(nrc):
            fcd = fc[j1]
            if j1 == 0:
                np1 = 0
            else:
                np1 = int(il[j1 - 1])
            np2 = int(il[j1] - 1)
            nrzt, mrks, mret = centra(np1, np2, fcd, b, c, epss, ksp, etp, blx, clx, epsc2[j1], a1[j1], a2[j1], nrzt, mrks, mret)
        
    else:
    # solução para o caso de flexo compressão
        for j1 in range(nrc):
            et01 = -c / b
            et12 = (-epsc2[j1] - c) / b
            fcd = fc[j1]
            if j1 == 0:
                np1 = 0
            else:
                np1 = int(il[j1 - 1])
            np2 = int(il[j1] - 1)
            for j2 in range(np1, np2):
                eps0 = epsp[j2]
                eps1 = epsp[j2 + 1]
                if eps0 == eps1:
                    continue
                if eps0 >= 0 and eps1 >= 0:
                    continue
                ks1i, et1i, ks2i, et2i, ks1ii, et1ii, ks2ii, et2ii = difer(j2, et01, et12, ksp, etp, eps0, eps1, epsc2[j1])
                nrzt, mrks, mret = regi(fcd, b, c, ks1i, et1i, ks2i, et2i, blx, clx, a1[j1], a2[j1], nrzt, mrks, mret)
                nrzt, mrks, mret = regii(fcd, ks1ii, et1ii, ks2ii, et2ii, nrzt, mrks, mret)

    nrz = nrzt
    mrx = mrks * ca - mret * sa
    mry = mrks * sa + mret * ca
    #print(nrz,mrx,mry, epss, epsi)

    return nrz, mrx, mry, epss, epsi

def difer(i, et01, et12, ksp, etp, eps0, eps1, epsc2):
    
    """
    Cálculo dos limites entre as regiões 0, I e II dentro de um segmento da poligonal
    """
    t01 = 0
    t12 = 0
    ks1i = 0
    et1i = 0
    ks2i = 0
    et2i = 0
    ks1ii = 0
    et1ii = 0
    ks2ii = 0
    et2ii = 0
    i2 = i + 1
    det = etp[i2] - etp[i]
    dksdet = (ksp[i2] - ksp[i]) / det
    dum1 = et01 - etp[i]
    dum2 = et12 - etp[i]
    ks01 = ksp[i] + dum1 * dksdet
    ks12 = ksp[i] + dum2 * dksdet
    det01 = dum1 / det
    det12 = dum2 / det
    if det01 > 0 and det01 < 1:
        t01 = 1
    if det12 > 0 and det12 < 1:
        t12 = 1
    if eps0 < eps1:
        t01 = -t01
        t12 = -t12
    if t01 == 0 and t12 == 0:
        if eps0 < 0:
            if eps0 > -epsc2:
                ks1i = ksp[i]
                et1i = etp[i]
                ks2i = ksp[i2]
                et2i = etp[i2]
            else:
                ks1ii = ksp[i]
                et1ii = etp[i]
                ks2ii = ksp[i2]
                et2ii = etp[i2]
    else:
        if t01 == 1:
            ks1i = ks01
            et1i = et01
            if t12 == 1:
                ks2i = ks12
                et2i = et12
                ks1ii = ks12
                et1ii = et12
                ks2ii = ksp[i2]
                et2ii = etp[i2]
            else:
                ks2i = ksp[i2]
                et2i = etp[i2]
        elif t01 == -1:
            ks2i = ks01
            et2i = et01
            if t12 == -1:
                ks1i = ks12
                et1i = et12
                ks2ii = ks12
                et2ii = et12
                ks1ii = ksp[i]
                et1ii = etp[i]
            else:
                ks1i = ksp[i]
                et1i = etp[i]
        else:
            if t12 == 1:
                ks1i = ksp[i]
                et1i = etp[i]
                ks2i = ks12
                et2i = et12
                ks1ii = ks12
                et1ii = et12
                ks2ii = ksp[i2]
                et2ii = etp[i2]
            else:
                ks1i = ks12
                et1i = et12
                ks2i = ksp[i2]
                et2i = etp[i2]
                ks1ii = ksp[i]
                et1ii = etp[i]
                ks2ii = ks12
                et2ii = et12
    return ks1i, et1i, ks2i, et2i, ks1ii, et1ii, ks2ii, et2ii
    


def centra(np1, np2, fcd, b, c, epss, ksp, etp, blx, clx, epsc2, a1, a2, nrzt, mrks, mret):
    """
    Solução para compressão centrada

    """
    if epss >= 0:
        return nrzt, mrks, mret
    if epss <= -epsc2:
        for j1 in range(np1, np2):
            j2 = j1 + 1
            regii(fcd, ksp[j1], etp[j1], ksp[j2], etp[j2], nrzt, mrks, mret)
    else:
        for j1 in range(np1, np2):
            j2 = j1 + 1
            regi(fcd, b, c, ksp[j1], etp[j1], ksp[j2], etp[j2], blx, clx, a1, a2, nrzt, mrks, mret)
    return nrzt, mrks, mret
    
       

def regi(fcd, b, c, ks1, et1, ks2, et2, blx, clx, a1, a2, nrzt, mrks, mret):
    """
    Solução para a integração dos esforços na região I
    """
    
    ble = 2.0 * a2 * b
    cle = 2.0 * a2 * c + a1
    d0 = c * a1 + a2 * c * c
    d1 = b * cle
    d2 = a2 * b * b
    e0 = cle * clx
    e1 = ble * clx + cle * blx
    e2 = ble * blx
    br = 0.85 * fcd
    dks = ks2 - ks1
    det = et2 - et1
    det1 = det / 2.0
    det2 = det * det
    det3 = det2 * det
    dks2 = dks * dks
    g00 = (ks1 + dks / 2.0) * det
    g01 = (ks1 * (et1 + det1) + dks * (et1 / 2.0 + det / 3.0)) * det
    g02 = (ks1 * (et1 * (det + et1) + det2 / 3.0) + dks * (et1 * (et1 / 2.0 + det / 1.5) + det2 / 4.0)) * det
    g03 = (ks1 * (et1 * (det2 + et1 * (1.5 * det + et1)) + det3 / 4.0) + dks * (et1 * (0.75 * det2 + et1 * (det + et1 / 2.0)) + det3 / 5.0)) * det
    g10 = (ks1 * (ks1 + dks) + dks2 / 3.0) * det1
    g11 = (ks1 * (ks1 * (et1 + det1) + dks * (et1 + det / 1.5)) + dks2 * (et1 / 3.0 + det / 4.0)) * det1
    g12 = (ks1 * (ks1 * (et1 * (et1 + det) + det2 / 3.0) + dks * (et1 * (et1 + det / 0.75) + det2 / 2.0))+ dks2 * (et1 * (et1 / 3.0 + det1) + det2 / 5.0)) * det1
    
    # cálculo dos esforços resistentes em relação aos eixos ksi e eta
    nrzt = nrzt + br * (d0 * g00 + d1 * g01 + d2 * g02)
    mrks = mrks + br * (d0 * g01 + d1 * g02 + d2 * g03)
    mret = mret - br * (d0 * g10 + d1 * g11 + d2 * g12)
    
    return nrzt, mrks, mret
    
    
def regii(fcd, ks1, et1, ks2, et2, nrzt, mrks, mret):
    """
    solução para a integração na região ii
    """
    det = et2 - et1
    dks = ks2 - ks1
    g00 = (ks1 + dks / 2.0) * det
    g01 = (ks1 * (et1 + det / 2.0) + dks * (et1 / 2.0 + det / 3.0)) * det
    g10 = (ks1 * (ks1 + dks) + dks * dks / 3.0) * det / 2.0
    fc = 0.85 * fcd

    # cálculo dos esforços resistentes em relação aos eixos ksi e eta
    nrzt = nrzt - fc * g00
    mrks = mrks - fc * g01
    mret = mret + fc * g10

    return nrzt, mrks, mret


def aco(e, epsb, fyd):
    """
    Calcula a tensão na barra de aço para uma dada deformação
    """
    eps2 = fyd / e
    if abs(epsb) <= eps2:
        sig = e * epsb
        et = e
    else:
        sig = abs(fyd) * (epsb / abs(epsb))
        et = 0
    return sig, et





#
# Início da solução do problema da flexão composta oblíqua
#
"""
Entrada de dados e chamada da rotina ajustl para o ajuste do equilíbrio entre esforços externos e internos

"""

arq = ""

while arq == "":
    print("\n>>>>>> qual o nome do arquivo de dados ?")
    arq = input().strip()

ir = 1
iw = 0
path_arq = "C://Users//Mauro//OneDrive//FCOpy//TAC2024//"
arq = path_arq + arq
with open(arq, 'r') as file:
    lines = file.readlines()

    # leitura da opção de cálculo: 1 = dimensionamento, 2 = verificação
    op = int(lines[ir - 1].strip())

    # coeficientes de segurança: gc = concreto, gs = aço
    gc, gs = map(float, lines[ir].split())

    # entrada da poligonal da seção de concreto: nv = número de vértices, xp e yp = coordenadas dos vértices
    nv = int(lines[ir + 1].strip())
    xp = np.zeros(nv)
    yp = np.zeros(nv)
    for i in range(nv):
        xp[i], yp[i] = map(float, lines[ir + 2 + i].split())

    # entrada no número de regiões da seção com concretos diferentes
    # nrc = número de regiões com concreto de diferente fck
    # fcd = vetor com o fck de cada região diferente
    # il = número do vértice em que a região de concreto termina
    nrc = int(lines[ir + 2 + nv].strip())
    fcd = np.zeros(nrc)
    il = np.zeros(nrc)
    for i in range(nrc):
        fcd[i], il[i] = map(float, lines[ir + 3 + nv + i ].split())
        fcd[i] /= gc  # converte fck em fcd dividindo por gamma_c
        il[i] = int(il[i]) # número do último vértice da região de concreto

    # leitura dos dados da armadura: nb = número de barras, e = módulo de elasticidade do aço, astotal = área de aço total
    astotal = 0.00
    if op == 1:
        nb, e = map(float, lines[ir + 3 + nv + nrc].split())
    else:
        nb, e, astotal = map(float, lines[ir + 3 + nv + nrc].split())
    nb = int(nb)
    xb = np.zeros(nb)
    yb = np.zeros(nb)
    fyd = np.zeros(nb)
    perc = np.zeros(nb)

    # leitura dos dados de cada barra: xb e yb = coordenadas, fyd = tensão de escoamento característica (fyk), perc = percentual de As,total
    for i in range(nb):
        xb[i], yb[i], fyd[i], perc[i] = map(float, lines[ir + 4 + nv + nrc + i].split())
        fyd[i] /= gs # converte fyk em fyd dividindo por gamma_s

    # esforços atuantes: na = esforço normal, max = momento fletor em torno do eixo x, may = momento fletor em torno do eixo y
    na, max, may = map(float, lines[ir + 4 + nv + nrc + nb].split())
    const = 210000.0 / e # constante para conversão de unidades para MPa
    
    # cálculo dos parâmetros da curva tensão deformação para cada região diferente do concreto
    epsc2 = np.zeros(nrc) # deformação correspondente ao final do trecho parabólico
    epscu = np.zeros(nrc) # deformação última de ruptura do concreto
    a1 = np.zeros(nrc) # coeficiente do termo do primeiro grau
    a2 = np.zeros(nrc) # coeficiente do termo do segundo grau

    # laços sobre as regiões com valores de fck diferentes
    for i in range(nrc):
        fck = fcd[i] * gc * const # calcula fck em MPa
        if fck <= 50:
            # parâmetros para concreto do grupo I: fck<= 50 MPa
            epsc2[i] = 0.002
            epscu[i] = 0.0035
            a1[i] = 1000.0
            a2[i] = 250000.0
        else:
            # parâmetros para concreto do grupo II: fck> 50 MPa
            epsc2[i] = 0.002 + 0.000085 * (fck - 50) ** 0.53
            epscu[i] = 0.0026 + 0.035 * ((90 - fck) / 100) ** 4
            nn = 1.4 + 23.4 * ((90 - fck) / 100) ** 4

            # cálculo dos coeficientes para a parábola do segundo grau equivalente à equação da norma com n diferente de 2
            ddx = epsc2[i] / 1000.0
            x = 0.0
            x2 = 0.0
            x3 = 0.0
            x4 = 0.0
            xy = 0.0
            x2y = 0.0
            while x <= epsc2[i]:
                y = 1.0 - (1.0 - x / epsc2[i]) ** nn
                x2 += x * x
                x3 += x * x * x
                x4 += x * x * x * x
                xy += x * y
                x2y += x * x * y
                x += ddx
            
            a1[i] = (x2y - x4 * xy / x3) / (x3 - x2 * x4 / x3)
            a2[i] = -(x2y - x3 * a1[i]) / x4

        # cálculo do coeficiente de fragilidade eta_c
        if fck <= 40:
            eta_c = 1.00
        else:
            eta_c = (40./fck)**(1./3.)
        fcd[i] = eta_c * fcd[i]
    
    # cálculo das propriedades geométricas da seção de concreto usando o Teorema de Green
        
    xmax = xp.max()
    xmin = xp.min()
    ymax = yp.max()
    ymin = yp.min()
    
    lx = xmax - xmin
    ly = ymax - ymin
    nv1 = nv - 1
    area = 0.0
    sx = 0.0
    sy = 0.0
    jx = 0.0
    jy = 0.0
    jxy = 0.0
    
    for i in range(nv1):
        dx = xp[i + 1] - xp[i]
        dy = yp[i + 1] - yp[i]
        area += (xp[i] + dx / 2.0) * dy
        sx += (xp[i] * (yp[i] + dy / 2.0) + dx * (yp[i] / 2.0 + dy / 3.0)) * dy
        sy += (xp[i] * (xp[i] + dx) + dx * dx / 3.0) * dy / 2.0
        jx += (xp[i] * (yp[i] * (dy + yp[i]) + dy * dy / 3.0) + dx * (yp[i] * (yp[i] / 2.0 + dy / 1.5) + dy * dy / 4.0)) * dy
        jy += (dx ** 3 / 4.0 + xp[i] * (dx * dx + xp[i] * (1.5 * dx + xp[i]))) * dy / 3.0
        jxy += (xp[i] * (xp[i] * (yp[i] + dy / 2.0) + dx * (yp[i] + dy / 1.5)) + dx * dx * (yp[i] / 3.0 + dy / 4.0)) * dy / 2.0
    
    # propriedades geométricas da seção transversal
    xg = sy / area # coordenada xg do centroide
    yg = sx / area # coordenada yg do centroide
    sox = sx - yg * area
    soy = sy - xg * area
    jox = jx - area * yg * yg   # momento de inércia em relação ao eixo x centroidal
    joy = jy - area * xg * xg   # momento de inércia em relação ao eixo y centroidal
    joxy = jxy - xg * yg * area # produto de inércia em relação sistema de coordenadas centroidal
    
    # converte as coordenadas dos vértices para o sistema centroidal de coordenadas
    xp = xp - xg
    yp = yp - yg
    
    
    # converte as coordenadas das barras de aço para o sistema centroidal de coordenadas
    xb = xb - xg
    yb = yb - yg

    # desenha a seção de concreto com as barras de aço
    # Criando a figura e os eixos
    fig, ax = plt.subplots()

    # Desenhando a seção de concreto
    largura = np.abs(xmax - xmin)
    altura = np.abs(ymax - ymin)
    ax.plot(xp,yp, '-')
    ax.plot(xb,yb, 'o',  color='black', markersize=5)


    # Ajustando a escala dos eixos para manter a proporção 1:2
    ax.set_aspect('equal')



    # Exibindo o gráfico
    plt.show()
    
    # Solução do sistema 3x3 não-linear de equações de equilíbrio
    """
    OP = 1 : Dimensionamento: ajusta as variáveis {x, alpha, As}, até que os esforços resistentes equilibrem os esforços externos
    OP = 2 : Verificação: calcula os esforços resistentes para um terno {x, alpha, As} inicial e imprime os resultados
    """
    
    # inicialização de variáveis
    nr, mrx, mry, lam, lammin, nrmin, mrxmin, mrymin = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
   
    dp = np.zeros(3)
    
    pi = math.pi
    pi2 = 1.5707963267948966192
    graus = 57.29577951308232

    # inicialização dos parâmetros da solução
    x = 0.00
    alph = 0.00
    tolmin = 1.0
    k = 0
    k0 = 0
    tol = 1000.00
    tole = 1e-8
    lam = 1.0 # fator de segurança inicial: lambda = 1.00
    b = 0.0 # zera a curvatura da seção
    c = 0.0 # zera a deformação normal no centroide da seção
    epss = 0.00
    epsi = 0.00


    x, alph, astotal = initial(jox, joy, joxy, lx, ly, na, max, may, fcd, fyd, astotal)
    
    if op == 1:
        # Solução do problema do dimensionamento

        x, alph, astotal = fsolve(dimen,(x, alph, astotal), args=(nv, xp, yp, nrc, fcd, il, nb, e , lam, xb, yb, fyd, perc, na, max, may,epsc2, epscu, a1, a2))
        nr, mrx, mry, epss, epsi = esfor(e, fyd, fcd, nv, xp, yp, nb, xb, yb, perc, x, alph, astotal, b, c, epss, epsi, nrc, il, epsc2, epscu, a1, a2)

    else:
        # Solução do problema da verificação
        x, alph, lam = fsolve(verif,(x, alph, lam), args=(nv, xp, yp, nrc, fcd, il, nb, e , astotal, xb, yb, fyd, perc, na, max, may,epsc2, epscu, a1, a2))
        nr, mrx, mry, epss, epsi = esfor(e, fyd, fcd, nv, xp, yp, nb, xb, yb, perc, x, alph, astotal, b, c, epss, epsi, nrc, il, epsc2, epscu, a1, a2)
    print("".ljust(78, "*") + "\n")
    print("**".ljust(76) + "**\n")

    if op == 2:
        print("**      verificacao de secao de concreto armado a solicitacoes normais      **\n")
    else:
        print("**     dimensionamento de secao de concreto armado a solicitacoes normais   **\n")

    print("**".ljust(76) + "**\n")
    print("".ljust(78, "*") + "\n")
    print("**                                           mx          my          n      **\n")
    print("**    esforcos atuantes de calculo:   {:12.4e}  {:12.4e}  {:12.4e}  **\n".format(max, may, na))
    print("**    esforcos resistentes de calculo:{:12.4e}  {:12.4e}  {:12.4e}  **\n".format(mrx, mry, nr))

    if op == 2:
        fs = 1.0 / lam
        print("**    fator de segurança:              fs = {:12.4e}                      **\n".format(fs))
    else:
        print("**             armadura:             area = {:12.4e}                      **\n".format(astotal))

    print("**    profundidade do linha neutra x:           {:12.4e}                  **\n".format(x))
    print("**    deformacao na fibra superior da secao:    {:12.4e}                  **\n".format(epss))
    print("**    deformacao na fibra inferior da secao:    {:12.4e}                  **\n".format(epsi))

    if abs(alph) > 2 * pi:
        alph = math.copysign(math.fmod(abs(alph), 2 * pi), alph)

    alpg = graus * alph
    print("**    inclinacao da linha neutra:              {:12.4e}                  **\n".format(alpg))
    print("".ljust(78, "*") + "\n")
    print("".ljust(78, "*") + "\n")
    print("\n")
    #print("".ljust(25) + "tecle <enter> para continuar\n\n")
    #car = input().strip()

    quit()


