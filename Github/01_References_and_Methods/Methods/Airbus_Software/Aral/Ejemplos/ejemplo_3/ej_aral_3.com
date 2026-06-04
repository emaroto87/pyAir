

2
1
ej_aral_3
     ! Unidades (separadas por una coma) (N/daN/kgf/lbf/klb , mm/cm/in)
N,mm
     ! Datos de los materiales. Dos posibilidades:
     !  MAT_ID,  MAT,FORM,TRAT,AREA,ESPE,BASE,DIRE : Fyp,Rp(%)
     !  MAT_ID [,MAT,FORM,TRAT,AREA,ESPE,BASE,DIRE]: E,Ec,G,µ,Fty,Fcy,Ftu,Fsu,e(%): Fyp,Rp(%)
MAT1
7075
Placas
T7651
0.
25.000
S
L
0
0
*
     ! Título de identificación del análisis (Máx. 80 caract.)
Ejemplo 3
     ! Datos de los larguerillos. Posibilidades:
     !  MATE_ID, L, C, T_Larg (0), CD o EX, T1, A, d, I, F_loc, F_crip, T_Unión (I0 P0 RSx o RTx)
     !  MATE_ID, L, C, T_Larg (1), CD o EX, B1, T1, B2, T2, T_Unión (P0 RSx o RTx)
     !  MATE_ID, L, C, T_Larg (2), B1, T1, B2, T2, T_Unión (P0 RSx o RTx)
     !  MATE_ID, L, C, T_Larg (3), CD o EX, B1, T1, B2, T2, B3, T3, T_Unión (P0 RSx o RTx)
     !  MATE_ID, L, C, T_Larg (4 o 5), B1, T1, B2, T2, B3, T3, T_Unión (P0 RSx o RTx)
     !  MATE_ID, L, C, T_Larg (6), B2, T2
     !  MATE_ID, L, C, T_Larg (7 u 8), B2, T2, B3, T3
MAT1
500
2
2
6
20
3
     ! Datos del revestimiento:
     !  MAT_ID, B1, T1, B2, T2, Tpp
MAT1
150
2.5
150
2
2.5
     ! Datos del "pad" y remachado:
     !  Bp, Tp [, Paso, Tipo_Rem (1 2 o 3)]
20
3
     ! Datos de cargas últimas:
     !  [Nom_Caso,] PC, SH1, SH2, p1, p2
CASO 1
-100000
0
0
0
0
CASO 2
-90000
20
10
0
0
CASO 3
-80000
40
30
0
0
CASO 4
-70000
60
50
0
0
CASO 5
-60000
80
70
0
0
CASO 6
-50000
100
90
0
0
CASO 7
-40000
120
110
0
0
CASO 8
-30000
140
130
0
0
CASO 9
-20000
160
150
0
0
CASO 10
-10000
180
170
0
0
CASO 11
0
180
170
0
0
*
