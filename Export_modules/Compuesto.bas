Attribute VB_Name = "Compuesto"
'Option Explicit

'Function M_Giro_Coord_Mat_to_Global(Alfa)

'Function M_Giro_Coord_Global_to_Mat(Alfa)

'Function M_Flex_Coord_Mat(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double)

'Function M_Rig_Coord_Mat(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double)

'Function M_Rig_Coord_Global(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double, Alfa As Double)

'Function M_Flex_Coord_Global(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double, Alfa As Double)

'M_K_Laminado(Eficiencia_range As Range, Alfa_range As Range, t_range As Range,_
'   E1_range As Range, E2_range As Range, G12_range As Range, Nu21_range As Range, Nu12_range As Range)

'Function R_Stress_Max(X, X_comp, Y, Y_comp, S, s1, s2, t12)

'Function R_TsaiHill(X, X_comp, Y, Y_comp, S, s1, s2, t12)

'Function Propiedades_Laminado(Matriz_Rigidez_Laminado As Range, t As Double)




'Funciones cambio de coordenadas

    Function M_Giro_Coord_Mat_to_Global(Alfa)
    'alfa en radianes
    'giro de coordenadas materiales a globales
        Dim c As Double
        Dim S As Double
        Dim M_Giro_Coord_Mat_to_Global_Temp() As Double
        ReDim M_Giro_Coord_Mat_to_Global_Temp(1 To 3, 1 To 3)
        
        c = Cos(Alfa)
        S = Sin(Alfa)
        
        M_Giro_Coord_Mat_to_Global_Temp(1, 1) = c ^ 2
        M_Giro_Coord_Mat_to_Global_Temp(1, 2) = S ^ 2
        M_Giro_Coord_Mat_to_Global_Temp(1, 3) = -2 * c * S
        M_Giro_Coord_Mat_to_Global_Temp(2, 1) = S ^ 2
        M_Giro_Coord_Mat_to_Global_Temp(2, 2) = c ^ 2
        M_Giro_Coord_Mat_to_Global_Temp(2, 3) = 2 * c * S
        M_Giro_Coord_Mat_to_Global_Temp(3, 1) = c * S
        M_Giro_Coord_Mat_to_Global_Temp(3, 2) = -c * S
        M_Giro_Coord_Mat_to_Global_Temp(3, 3) = c ^ 2 - S ^ 2
        
        M_Giro_Coord_Mat_to_Global = M_Giro_Coord_Mat_to_Global_Temp
    
    End Function
    
    Function M_Giro_Coord_Global_to_Mat(Alfa)
    'alfa en radianes
    'giro de coordenadas materiales a globales
        
        M_Giro_Coord_Global_to_Mat = Application.WorksheetFunction.MInverse(M_Giro_Coord_Mat_to_Global(Alfa))
    
    End Function



'Matrices Lámina

    Function M_Flex_Coord_Mat(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double)
        
        Dim M_Flex_Coord_Mat_Temp() As Double
        
        ReDim M_Flex_Coord_Mat_Temp(1 To 3, 1 To 3)
        
        M_Flex_Coord_Mat_Temp(1, 1) = 1 / E1
        M_Flex_Coord_Mat_Temp(1, 2) = -Nu12 / E2
        M_Flex_Coord_Mat_Temp(1, 3) = 0
        M_Flex_Coord_Mat_Temp(2, 1) = -Nu12 / E2
        M_Flex_Coord_Mat_Temp(2, 2) = 1 / E2
        M_Flex_Coord_Mat_Temp(2, 3) = 0
        M_Flex_Coord_Mat_Temp(3, 1) = 0
        M_Flex_Coord_Mat_Temp(3, 2) = 0
        M_Flex_Coord_Mat_Temp(3, 3) = 1 / G12
        
        M_Flex_Coord_Mat = M_Flex_Coord_Mat_Temp
    
    End Function
    
    
    Function M_Rig_Coord_Mat(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double)
        Dim M_Flex_Coord_Mat_Temp() As Double
        ReDim M_Flex_Coord_Mat_Temp(1 To 3, 1 To 3)
    
        M_Flex_Coord_Mat_Temp = M_Flex_Coord_Mat(E1, E2, G12, Nu21, Nu12)
        M_Rig_Coord_Mat = Application.WorksheetFunction.MInverse(M_Flex_Coord_Mat_Temp)
    
    End Function
    
    
    Function M_Rig_Coord_Global(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double, Alfa As Double)
        Dim M_Rig_Coord_Mat_Temp As Variant
        Dim t As Variant
        Dim T_trasp() As Variant
        Dim M_Rig_Coord_Global_Temp As Variant
        
        t = M_Giro_Coord_Mat_to_Global(Alfa)
        
        M_Rig_Coord_Mat_Temp = M_Rig_Coord_Mat(E1, E2, G12, Nu21, Nu12)
    
        With Application.WorksheetFunction
            T_trasp = .Transpose(t)
            M_Rig_Coord_Global_Temp = .MMult(t, M_Rig_Coord_Mat_Temp)
            M_Rig_Coord_Global_Temp = .MMult(M_Rig_Coord_Global_Temp, T_trasp)
        End With
        
        M_Rig_Coord_Global = M_Rig_Coord_Global_Temp
    
    End Function
    
    Function M_Flex_Coord_Global(E1 As Double, E2 As Double, G12 As Double, Nu21 As Double, Nu12 As Double, Alfa As Double)
        Dim M_Rig_Coord_Global_Temp() As Variant
    
        M_Rig_Coord_Global_Temp = M_Rig_Coord_Global(E1, E2, G12, Nu21, Nu12, Alfa)
        M_Flex_Coord_Global = Application.WorksheetFunction.MInverse(M_Rig_Coord_Global_Temp)
    End Function



'Matrices Laminado

    Function M_K_Laminado(Eficiencia_range As Range, Alfa_range As Range, t_range As Range, _
        E1_range As Range, E2_range As Range, G12_range As Range, Nu21_range As Range, Nu12_range As Range)
        
        Dim Eficiencia() As Double
        Dim Alfa() As Double
        Dim t() As Double
        Dim E1() As Double
        Dim E2() As Double
        Dim G12() As Double
        Dim Nu21() As Double
        Dim Nu12() As Double
        
        Dim i As Integer
        Dim j As Integer
        Dim k As Integer
        
        Dim Dimension As Integer 'dimensión de los vectores de entrada
        Dim n_Laminas As Integer
        Dim t_Total As Double
        
        Dim z0_Fibra() As Double 'coordenada inicial de la lámina i
        Dim z1_Fibra() As Double ' coordenada final de la lámina i
        
    
        
    'número de láminas
        Dimension = UBound(Alfa_range(), 1)
        For i = 1 To Dimension
            If Alfa_range(i).Value = "" Then Exit For
        Next
        n_Laminas = i - 1
        
    'paso los valores del los rangos a matrices de valores para no tener que estar poniendo constantemente vect(i).value
        ReDim Eficiencia(1 To n_Laminas)
        ReDim Alfa(1 To n_Laminas)
        ReDim t(1 To n_Laminas)
        ReDim E1(1 To n_Laminas)
        ReDim E2(1 To n_Laminas)
        ReDim G12(1 To n_Laminas)
        ReDim Nu21(1 To n_Laminas)
        ReDim Nu12(1 To n_Laminas)
    
        For i = 1 To n_Laminas
            Eficiencia(i) = Eficiencia_range(i).Value
            Alfa(i) = Alfa_range(i).Value
            t(i) = t_range(i).Value
            E1(i) = E1_range(i).Value
            E2(i) = E2_range(i).Value
            G12(i) = G12_range(i).Value
            Nu21(i) = Nu21_range(i).Value
            Nu12(i) = Nu12_range(i).Value
        Next
    
    'espesor total
        For i = 1 To n_Laminas
            t_Total = t(i) + t_Total
        Next
    
    'coordenadas Z
        ReDim z0_Fibra(1 To n_Laminas)
        ReDim z1_Fibra(1 To n_Laminas)
    
        For i = 1 To n_Laminas
            If i = 1 Then
                z0_Fibra(1) = t_Total / 2
            Else
                z0_Fibra(i) = z1_Fibra(i - 1)
            End If
            z1_Fibra(i) = z0_Fibra(i) - t(i)
        Next
    
    'matrices A, B, C
        Dim M_Rigidez_Lamina As Variant
        Dim M_A(1 To 3, 1 To 3) As Double
        Dim M_B(1 To 3, 1 To 3) As Double
        Dim M_D(1 To 3, 1 To 3) As Double
'        ReDim M_Rigidez_Lamina(1 To 3, 1 To 3)
        
        For i = 1 To n_Laminas
            If Eficiencia(i) = 0 Then GoTo Siguiente_Lamina
            M_Rigidez_Lamina = M_Rig_Coord_Global(E1(i), E2(i), G12(i), Nu21(i), Nu12(i), Alfa(i))
            For j = 1 To 3
                For k = 1 To 3
                    M_A(j, k) = M_A(j, k) + (z0_Fibra(i) - z1_Fibra(i)) * M_Rigidez_Lamina(j, k)
                    M_B(j, k) = M_B(j, k) + (1 / 2) * (z0_Fibra(i) ^ 2 - z1_Fibra(i) ^ 2) * M_Rigidez_Lamina(j, k)
                    M_D(j, k) = M_D(j, k) + (1 / 3) * (z0_Fibra(i) ^ 3 - z1_Fibra(i) ^ 3) * M_Rigidez_Lamina(j, k)
                Next
            Next
Siguiente_Lamina:
        Next
    
    'matriz K del laminado
        Dim M_K_Laminado_temp(1 To 6, 1 To 6)
        For i = 1 To 3
            For j = 1 To 3
                M_K_Laminado_temp(i, j) = M_A(i, j)
                M_K_Laminado_temp(i, j + 3) = M_B(i, j)
                M_K_Laminado_temp(i + 3, j) = M_B(i, j)
                M_K_Laminado_temp(i + 3, j + 3) = M_D(i, j)
            Next
        Next
    
        M_K_Laminado = M_K_Laminado_temp
    
    End Function


'Rotura
    Function R_Stress_Max(X, X_comp, Y, Y_comp, S, s1, s2, t12)
        Dim R1 As Double
        Dim R2 As Double
        Dim R12 As Double
                
        If s1 < 0 Then
            R1 = -s1 / X_comp
        Else
            R1 = s1 / X
            
        End If
        
        If s2 < 0 Then
            R2 = -s2 / Y_comp
        Else
            R2 = s2 / Y
        End If
        
        R12 = Abs(t12) / S
        
        R_Stress_Max = Application.WorksheetFunction.MAX(R1, R2, R12)
    End Function
    
    
    Function R_TsaiHill(X, X_comp, Y, Y_comp, S, s1, s2, t12)
        R_TsaiHill = (s1 / X) ^ 2 + (s2 / Y) ^ 2 + (t12 / S) ^ 2 - s1 * s2 / X ^ 2
    End Function


'Propiedades elásticas del laminado
    Function Propiedades_Laminado(Matriz_Rigidez_Laminado As Range, t As Double)
    Dim Matriz_Flexibilidad_Laminado As Variant
    Dim Exx As Double
    Dim Eyy As Double
    Dim Gxy As Double
    Dim nu_xy As Double
    Dim nu_yx As Double
    Dim A11_Inv As Double
    Dim A22_Inv As Double
    Dim A33_Inv As Double
    Dim A12_Inv As Double
    Dim Propiedades_Laminado_temp(1 To 5) As Double
    
    Matriz_Flexibilidad_Laminado = Application.WorksheetFunction.MInverse(Matriz_Rigidez_Laminado)
    
    A11_Inv = Matriz_Flexibilidad_Laminado(1, 1)
    A22_Inv = Matriz_Flexibilidad_Laminado(2, 2)
    A33_Inv = Matriz_Flexibilidad_Laminado(3, 3)
    A12_Inv = Matriz_Flexibilidad_Laminado(1, 2)
    
    Exx = 1 / (t * A11_Inv)
    Eyy = 1 / (t * A22_Inv)
    Gxy = 1 / (t * A33_Inv)
    nu_xy = -A12_Inv / A11_Inv
    nu_yx = -A12_Inv / A22_Inv
       
    Propiedades_Laminado_temp(1) = Exx
    Propiedades_Laminado_temp(2) = Eyy
    Propiedades_Laminado_temp(3) = Gxy
    Propiedades_Laminado_temp(4) = nu_xy
    Propiedades_Laminado_temp(5) = nu_yx

    Propiedades_Laminado = Propiedades_Laminado_temp
    
    End Function
