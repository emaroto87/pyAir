Attribute VB_Name = "FASTENERS"
Function interaction(Shear_Applied As Double, Shear_Strength As Double, Shear_Coeff As Double, Axial_Applied As Double, Axial_Strength As Double, Axial_Coeff As Double) As Double
    'Función que calcula el RF de interacción en un fástener sometido al mismo tiempo a cortadura y tracción
    Dim RF As Double
    RF = bolzano("((x*" & CStr(Shear_Applied) & ")/" & CStr(Shear_Strength) & ")^" & CStr(Shear_Coeff) & "+" & "((x*" & CStr(Axial_Applied) & ")/" & CStr(Axial_Strength) & ")^" & CStr(Axial_Coeff) & "-1", 0, 1000000, 0.000001)(0)
    interaction = RF
End Function

Function K_Tronco_Cono(Young As Double, Alfa As Double, dint As Double, Dext_ini As Double, ti As Double, tf As Double) As Double
    ' Función que devuelve la rigidez de un tronco de cono
    Dim pii As Double
    pii = WorksheetFunction.pi()
    Dim tangente As Double
    tangente = Tan(Alfa * pii / 180)
    K_Tronco_Cono = (pii * Young * dint * tangente) / (Log((2 * tf * tangente - dint + Dext_ini) / (2 * tf * tangente + dint + Dext_ini)) - Log((2 * ti * tangente - dint + Dext_ini) / (2 * ti * tangente + dint + Dext_ini)))
End Function
Function K_Tronco_Cono1(Young As Double, Alfa As Double, dint As Double, Dext_ini As Double, t As Double) As Double
    ' Función que devuelve la rigidez de un tronco de cono
    Dim pii As Double
    pii = WorksheetFunction.pi()
    Dim tangente As Double
    tangente = Tan(Alfa * pii / 180)
    K_Tronco_Cono1 = (pii * Young * dint * tangente) / Log(((2 * t * tangente - dint + Dext_ini) * (Dext_ini + dint)) / ((2 * t * tangente + dint + Dext_ini) * (Dext_ini - dint)))
End Function

Function K_Cilindro_Hueco(Young As Double, dint As Double, Dext As Double, ti As Double, tf As Double) As Double
    'Función que devuelve la rigidez de un cilindro hueco
    Dim pii As Double
    pii = WorksheetFunction.pi()
    K_Cilindro_Hueco = Young * pii * ((0.5 * Dext) ^ 2 - (0.5 * dint) ^ 2) / (tf - ti)
End Function

Function K_Bolt(young_bolt As Double, Shank_diameter As Double, Bolt_Head_Height As Double, Nut_Height As Double, Stack_Length As Double) As Double
    ' Función que devuelve la rigidez de un fástener
    Dim pii As Double
    pii = WorksheetFunction.pi()
    K_Bolt = young_bolt * pii * (0.5 * Shank_diameter) ^ 2 / (Bolt_Head_Height / 3 + Nut_Height / 3 + Stack_Length)
End Function

Function Washer_Diameter(Head_Ext_Diameter As Double, Washer_Ext_Diameter As Double, Washer_Thickness, Alfa As Double) As Double
    ' Función que devuelve el diámetro efectivo de una arandela
    Dim pii As Double
    pii = WorksheetFunction.pi()
    Dim tangente As Double
    tangente = Tan(Alfa * pii / 180)

    If Washer_Ext_Diameter = 0 Or Washer_Thickness = 0 Then
        Washer_Diameter = Head_Ext_Diameter
    Else
        Dim delta As Double
        delta = Washer_Thickness * tangente
        
        If (Head_Ext_Diameter + 2 * delta) > Washer_Ext_Diameter Then
            Washer_Diameter = Washer_Ext_Diameter
        Else
            Washer_Diameter = Head_Ext_Diameter + 2 * delta
        End If
    End If
End Function

Function Plano_Medio_Stack(Bolt_Head_Ext_Dia As Double, Bolt_Washer_Ext_Dia As Double, Bolt_Washer_Thickness As Double, Alfa As Double, Nut_Head_Ext_Dia As Double, Nut_Washer_Ext_Dia As Double, Nut_Washer_Thickness As Double, Stack As Range)
    ' Función que devuelve la posición del plano medio de un apliamiento de placas unidas por un fástener
    Dim Washer_Effective_Diameter_Bolt As Double
    Washer_Effective_Diameter_Bolt = Washer_Diameter(Bolt_Head_Ext_Dia, Bolt_Washer_Ext_Dia, Bolt_Washer_Thickness, Alfa)
    
    Dim Washer_Effective_Diameter_Nut As Double
    Washer_Effective_Diameter_Nut = Washer_Diameter(Nut_Head_Ext_Dia, Nut_Washer_Ext_Dia, Nut_Washer_Thickness, Alfa)
    
    Dim pii As Double
    pii = WorksheetFunction.pi()
   
    Dim Stack_Thickness As Double
    Stack_Thickness = WorksheetFunction.Sum(Stack)
    
    Dim resultados(2) As Double
    resultados(0) = (1 / 4) * Tan((Alfa - 90) * pii / 180) * (Washer_Effective_Diameter_Nut - Washer_Effective_Diameter_Bolt) - Bolt_Washer_Thickness - (1 / 2) * Stack_Thickness
    resultados(1) = (1 / 4) * (Washer_Effective_Diameter_Bolt + Washer_Effective_Diameter_Nut) - (1 / 2) * Stack_Thickness / Tan((Alfa - 90) * pii / 180)
    resultados(2) = Stack_Thickness
    Plano_Medio_Stack = resultados
End Function

Function K_Stack(Bolt_Head_Ext_Dia As Double, Bolt_Washer_Ext_Dia As Double, Bolt_Washer_Thickness As Double, Bolt_Washer_Young As Double, Alfa As Double, Nut_Head_Ext_Dia As Double, Nut_Washer_Ext_Dia As Double, Nut_Washer_Thickness As Double, Nut_Washer_Young As Double, Stack_Thickness As Range, Stack_Efficiency As Range, Stack_Young As Range, Shank_Dia As Double, Bore_Dia As Double, Part_Max_With_Relative_To_Bore_Dia As Double, selector As String)
    ' Función que devuelve la rigidez de un apilamiento de placas unidas por un fástener (arandelas incluidas)
    Dim n As Double
    Dim m As Double
    Dim pii As Double
    Dim testigo As Integer
    testigo = 0
    Dim temporal1 As Double
    Dim temporal2 As Double
    pii = WorksheetFunction.pi()

    ' Cálculo de la posición del plano medio con respecto a la base de la cabeza del tornillo
    Dim Mid_Plane_Position As Double
    Mid_Plane_Position = Plano_Medio_Stack(Bolt_Head_Ext_Dia, Bolt_Washer_Ext_Dia, Bolt_Washer_Thickness, Alfa, Nut_Head_Ext_Dia, Nut_Washer_Ext_Dia, Nut_Washer_Thickness, Stack_Thickness)(0)
    
    ' Cálculo del máximo ancho de la distribución de tensión
    Dim Part_Max_Width As Double
    Part_Max_Width = Plano_Medio_Stack(Bolt_Head_Ext_Dia, Bolt_Washer_Ext_Dia, Bolt_Washer_Thickness, Alfa, Nut_Head_Ext_Dia, Nut_Washer_Ext_Dia, Nut_Washer_Thickness, Stack_Thickness)(1)
    
    ' Cálculo del máximo ancho de la distribución de tensión
    Dim Total_Stack_Thickness As Double
    Total_Stack_Thickness = Plano_Medio_Stack(Bolt_Head_Ext_Dia, Bolt_Washer_Ext_Dia, Bolt_Washer_Thickness, Alfa, Nut_Head_Ext_Dia, Nut_Washer_Ext_Dia, Nut_Washer_Thickness, Stack_Thickness)(2)
    
    
    ' Cálculo del máximo ancho de la distribución de tensión admisible. Se calcula a partir de un valor (introducido por el usuario) de ancho máximo admisible con respecto al agujero
    ' Si el máximo ancho de la distribución supera este valor admisible se trunca la zona cónica de la distribución y se convierte en una zona cilíndrica
    ' Si no se quiere tener en cuenta este efecto bastará con introducir un valor "infinito "de "Part_Max_With_Relative_To_Bore_Dia" (NORMALMENTE 2.5)
    Dim Part_Max_Width_Allowable As Double
    Part_Max_Width_Allowable = (Part_Max_With_Relative_To_Bore_Dia * Bore_Dia) / 2
    
    ' Se cuenta el número de placas que hay en el apilamiento
    Dim Stack_Number As Double
    Stack_Number = WorksheetFunction.CountIf(Stack_Efficiency, "1")
    
    ' Se cuenta la longitud del rango para poder recorrerlo
    Dim Range_Length As Double
    Range_Length = maxrange(Stack_Efficiency)
    
    ' Si el máximo de la distribución supera el admisible se calculan los dos puntos de intersección con las rectas que forman la zona cónica
    'If Part_Max_Width > Part_Max_Width_Allowable Then
        Dim Y_Bolt As Double
        Dim Y_Nut As Double
        Dim WEDB As Double ' Washer_Efective_Diameter_Bolt
        Dim WEDN As Double ' Washer_Efective_Diameter_Nut
        WEDB = Washer_Diameter(Bolt_Head_Ext_Dia, Bolt_Washer_Ext_Dia, Bolt_Washer_Thickness, Alfa)
        WEDN = Washer_Diameter(Nut_Head_Ext_Dia, Nut_Washer_Ext_Dia, Nut_Washer_Thickness, Alfa)
        Y_Bolt = Tan((Alfa - 90) * pii / 180) * (Part_Max_Width_Allowable - (1 / 2) * WEDB) - Bolt_Washer_Thickness ' Punto de inicio de la zona cilíndrica
        ST = WorksheetFunction.SumProduct(Stack_Efficiency, Stack_Thickness) ' Espesor total del apilamiento (sin arandelas) teniendo en cuenta la eficiencia. Se necesita para la fórmula de Y_Nut (ver teoría)
        Y_Nut = Tan((90 - Alfa) * pii / 180) * (Part_Max_Width_Allowable - (1 / 2) * WEDN) - Bolt_Washer_Thickness - ST ' Punto de fin de la zona cilíndrica
    If Part_Max_Width > Part_Max_Width_Allowable Then
        testigo = 1
        ' Se crea un vector que almacene la posición de el cambio de placa así como de los puntos que definen la zona cilíndrica
        ' Al redimensionar al número de placas + 1 + la posición cero tenemos (n+2) posiciones que necesitamos
        ReDim vec1(Stack_Number + 1) As Double
        ' Se rellena el vector
        m = 0
        For n = 1 To Range_Length
            If Stack_Efficiency(n, 1) = 1 Then
                If m = 0 Then
                    vec1(m) = Stack_Thickness(n, 1)
                Else
                    vec1(m) = vec1(m - 1) + Stack_Thickness(n, 1)
                End If
                m = m + 1
            End If
        Next n
        ' Por último, se colocan los puntos que forman la zona cilíndrica
        vec1(m) = Abs(Y_Bolt) - Bolt_Washer_Thickness
        vec1(m + 1) = Abs(Y_Nut) - Bolt_Washer_Thickness
        ' Podría ocurrir que la posición de inicio y fin de la zona cilíndrica coincida con un cambio de placas. Para arreglarlo se eliminan duplicados
        ' Se hace paso a paso creando distintos vectores "vecx" para facilitar la comprensión y trazabilidad de la función
        ReDim vec2(Stack_Number + 1) As Double
        vec2 = duplicados_vect(vec1)
        ' Si lo anterior hubiera ocurrido aparecerán uno o dos ceros en las últimas posiciónes del vector. En ese caso, se redimensiona el vector con una o dos posiciónes menos (ya de paso, para aprovechar el "if" se crea el siguiente vector)
        ReDim vec3(Stack_Number + 1) As Double
        vec3 = vec2
        

        If vec3(Stack_Number) = 0 Then
            ReDim Preserve vec3(Stack_Number - 1)
            ReDim vec4(Stack_Number - 1) As Double
        ElseIf vec3(Stack_Number + 1) = 0 Then
            ReDim Preserve vec3(Stack_Number)
            ReDim vec4(Stack_Number) As Double
        Else
            ReDim vec4(Stack_Number + 1) As Double
        End If
        ' Se ordena el vector de menor a mayor
        vec4 = ordena_vect(vec3)

        
    Else
        ' Se crea un vector que almacene la posición de el cambio de placa así como del punto de intersección de las rectas que forman las zonas cónicas
        ' Al redimensionar al número de placas + la posición cero tenemos (n+1) posiciones que necesitamos
        ReDim vec1(Stack_Number) As Double
        ' Se rellena el vector
        m = 0
        For n = 1 To Range_Length
            If Stack_Efficiency(n, 1) = 1 Then
                If m = 0 Then
                    vec1(m) = Stack_Thickness(n, 1)
                Else
                    vec1(m) = vec1(m - 1) + Stack_Thickness(n, 1)
                End If
                m = m + 1
            End If
        Next n
        ' Por último, se coloca el punto de intersección de las rectas que forman las zonas cónicas
        vec1(m) = Abs(Mid_Plane_Position) - Bolt_Washer_Thickness
        ' Podría ocurrir que la posición de la intersección de las rectas que forman las zonas cónicas zonas cónicas coincida con un cambio de placas. Para arreglarlo se eliminan duplicados
        ' Se hace paso a paso creando distintos vectores "vecx" para facilitar la comprensión y trazabilidad de la función
        ReDim vec2(Stack_Number) As Double
        vec2 = duplicados_vect(vec1)
        ' Si lo anterior hubiera ocurrido aparecerá un cero en la última posición del vector. En ese caso, se redimensiona el vector con una posición menos (ya de paso, para aprovechar el "if" se crea el siguiente vector)
        ReDim vec3(Stack_Number) As Double
        vec3 = vec2
        If vec3(Stack_Number) = 0 Then
            ReDim Preserve vec3(Stack_Number - 1)
            ReDim vec4(Stack_Number - 1) As Double
        Else
            ReDim vec4(Stack_Number) As Double
        End If
        ' Se ordena el vector de menor a mayor
        vec4 = ordena_vect(vec3)
    End If
    ' Independientemente del camino anterior, ahora existe un vector "vec4" que contiene ordenados (y sin huecos vacíos) todos los puntos de cambio de integración
    ' Se crea la matriz que contendrá toda la información necesaria para obtener las rigideces de cada tramo
    
    ' Se mide el número final de tramos del stack
    Dim tramos As Integer
    tramos = UBound(vec4(), 1) + 1
    
    ' Se crea la tabla de compilación:
    ReDim tabla(1 To tramos, 1 To 10) As Double
    
    ' Se rellena la columna 1 de la tabla: Número de orden del tramo
    For n = 1 To tramos
        tabla(n, 1) = n
    Next n
    
    ' Se rellena la columna 2 de la tabla: Posición inicial del tramo
    tabla(1, 2) = 0
    For n = 2 To tramos
        tabla(n, 2) = vec4(n - 2)
    Next n
    
    ' Se rellena la columna 3 de la tabla: Posición final del tramo
    For n = 1 To tramos
        tabla(n, 3) = vec4(n - 1)
    Next n
    
    ' Se rellena la columna 4 de la tabla: Posición media del tramo
    For n = 1 To tramos
        tabla(n, 4) = (tabla(n, 2) + tabla(n, 3)) / 2
    Next n
    
    ' Se rellena la columna 5 de la tabla: Diámetro del Bore
    For n = 1 To tramos
        tabla(n, 5) = Bore_Dia
    Next n
    
    ' Se rellena la columna 6 de la tabla: Ángulo
    For n = 1 To tramos
        tabla(n, 6) = Alfa
    Next n

    ' Se rellena la columna 7 de la tabla: Radio inicial del tramo
    For n = 1 To tramos
    
        temporal1 = WEDB / 2 + (Bolt_Washer_Thickness - tabla(n, 2) - Bolt_Washer_Thickness) / (Tan((tabla(n, 6) - 90) * pii / 180))
        temporal2 = WEDN / 2 + (Bolt_Washer_Thickness - tabla(n, 2) - Bolt_Washer_Thickness + ST) / (Tan((90 - tabla(n, 6)) * pii / 180))
        
        If testigo = 0 And temporal1 <= Part_Max_Width Then
            tabla(n, 7) = temporal1
        End If
        
        If testigo = 0 And temporal1 > Part_Max_Width Then
            tabla(n, 7) = temporal2
        End If
        
        If testigo = 1 And temporal1 <= Part_Max_Width_Allowable And (-tabla(n, 2) - Bolt_Washer_Thickness) > Mid_Plane_Position Then
            tabla(n, 7) = temporal1
        End If
        
        If testigo = 1 And temporal1 > Part_Max_Width_Allowable And temporal2 > Part_Max_Width_Allowable Then
            tabla(n, 7) = Part_Max_Width_Allowable
        End If
        
        If testigo = 1 And temporal2 <= Part_Max_Width_Allowable And (-tabla(n, 2) - Bolt_Washer_Thickness) < Mid_Plane_Position Then
            tabla(n, 7) = temporal2
        End If
        
    Next n

    ' Se rellena la columna 8 de la tabla: Radio final del tramo
    For n = 1 To tramos
    
        temporal1 = WEDB / 2 + (Bolt_Washer_Thickness - tabla(n, 3) - Bolt_Washer_Thickness) / (Tan((tabla(n, 6) - 90) * pii / 180))
        temporal2 = WEDN / 2 + (Bolt_Washer_Thickness - tabla(n, 3) - Bolt_Washer_Thickness + ST) / (Tan((90 - tabla(n, 6)) * pii / 180))
        
        If testigo = 0 And temporal1 <= Part_Max_Width Then
            tabla(n, 8) = temporal1
        End If
        
        If testigo = 0 And temporal1 > Part_Max_Width Then
            tabla(n, 8) = temporal2
        End If
        
        If testigo = 1 And temporal1 <= Part_Max_Width_Allowable And (-tabla(n, 3) - Bolt_Washer_Thickness) > Mid_Plane_Position Then
            tabla(n, 8) = temporal1
        End If
        
        If testigo = 1 And temporal1 > Part_Max_Width_Allowable And temporal2 > Part_Max_Width_Allowable Then
            tabla(n, 8) = Part_Max_Width_Allowable
        End If
        
        If testigo = 1 And temporal2 <= Part_Max_Width_Allowable And (-tabla(n, 3) - Bolt_Washer_Thickness) < Mid_Plane_Position Then
            tabla(n, 8) = temporal2
        End If
        
    Next n

    ' Se rellena la columna 9 de la tabla: Modulo de elasticidad "E"
    Dim acumulado As Double
    For n = 1 To tramos
        m = 1
        'acumulado = Stack_Thickness(1, 1)
        acumulado = 0
        Do While tabla(n, 3) > acumulado
            If Stack_Efficiency(m, 1) = 1 Then
                acumulado = acumulado + Stack_Thickness(m, 1)
            End If
                m = m + 1
        Loop
       tabla(n, 9) = Stack_Young(m - 1, 1)
        'tabla(n, 9) = 71000 'test
    Next n
    
    ' Se rellena la columna 10 de la tabla: Cálculo de la rigidez del tramo
    For n = 1 To tramos
    
        If tabla(n, 8) - tabla(n, 7) > 0.0000001 Then  'Actualización 25-01-2019 --> Pongo también el error aquí en lugar de una comparación directa por si los logaritmos dieran problemas. Ver explicación más abajo
            tabla(n, 10) = K_Tronco_Cono(tabla(n, 9), tabla(n, 6), tabla(n, 5), 2 * tabla(n, 7), tabla(n, 2), tabla(n, 3))
            'tabla(n, 10) = K_Tronco_Cono1(tabla(n, 9), tabla(n, 6), tabla(n, 5), 2 * tabla(n, 7), tabla(n, 3) - tabla(n, 2))
        End If
        
        If tabla(n, 7) - tabla(n, 8) > 0.0000001 Then  'Actualización 25-01-2019 --> Pongo también el error aquí en lugar de una comparación directa por si los logaritmos dieran problemas. Ver explicación más abajo
            tabla(n, 10) = K_Tronco_Cono(tabla(n, 9), tabla(n, 6), tabla(n, 5), 2 * tabla(n, 8), Total_Stack_Thickness - tabla(n, 3), Total_Stack_Thickness - tabla(n, 2))
            'tabla(n, 10) = K_Tronco_Cono1(tabla(n, 9), tabla(n, 6), tabla(n, 5), 2 * tabla(n, 8), tabla(n, 3) - tabla(n, 2))
        End If

        'Dependiendo de los datos de entrada, por culpa de la función "pi()" y de la función "Tan()", ocurre que aparece un error de 10^-16 que hace que el error no se cumpla y se calcule
        'como que tabla(n, 8)<tabla(n, 7) o que tabla(n, 8)<tabla(n, 7) --> En estos casos, aunque el valor de tabla(n, 10) se calcule mal en uno de los dos "If" anteriores aquí se sobreescribirá
        'Por este motivo, no se comprueba la igualdad, sino que cumplen con una condición de error. Este valor de error está fijado en 10^-7 que es más que suficiente para que no interfiera
        'con espesores razonables de placa
        Dim check As Double
        check = tabla(n, 8) - tabla(n, 7)
        If Abs(tabla(n, 8) - tabla(n, 7)) < 0.0000001 Then
            tabla(n, 10) = K_Cilindro_Hueco(tabla(n, 9), tabla(n, 5), 2 * tabla(n, 7), tabla(n, 2), tabla(n, 3))
        End If
        
    Next n
    
    
    '-----------------------------------------
    'Calculo de la rigidez del stack de piezas sin arandelas
    Dim K_Stack_without_washers As Double
    K_Stack_without_washers = 0
    For n = 1 To tramos
        K_Stack_without_washers = K_Stack_without_washers + 1 / tabla(n, 10)
    Next n
    K_Stack_without_washers = 1 / K_Stack_without_washers
    '-----------------------------------------
    
    
    
    'Variables auxiliares para que el cálculo de la rigidez de las arandelas quede más desglosado y trazable
    Dim temp1 As Double
    Dim temp2 As Double
    Dim temp3 As Double
    Dim temp4 As Double
    
    
    '-----------------------------------------
    'Calculo de la rigidez de la arandela de la cabeza del fástener
    Dim K_Washer_Bolt As Double
    If Bolt_Washer_Thickness = 0 Or Bolt_Washer_Ext_Dia = 0 Then
        K_Washer_Bolt = 10 ^ 100
    Else
        If Bolt_Washer_Ext_Dia > WEDB Then
            K_Washer_Bolt = K_Tronco_Cono1(Bolt_Washer_Young, Alfa, Bore_Dia, Bolt_Head_Ext_Dia, Bolt_Washer_Thickness)
        ElseIf Bolt_Head_Ext_Dia >= WEDB Then
            K_Washer_Bolt = K_Cilindro_Hueco(Bolt_Washer_Young, Bore_Dia, WEDB, 0, Bolt_Washer_Thickness)
        Else
            'Cálculo de la parte cónica
            temp1 = K_Tronco_Cono1(Bolt_Washer_Young, Alfa, Bore_Dia, Bolt_Head_Ext_Dia, (0.5 * (WEDB - Bolt_Head_Ext_Dia)) / (Tan((Alfa) * pii / 180)))
            'Cálculo de la parte cilindrica
            temp2 = K_Cilindro_Hueco(Bolt_Washer_Young, Bore_Dia, WEDB, 0, Bolt_Washer_Thickness - (0.5 * (WEDB - Bolt_Head_Ext_Dia)) / (Tan((Alfa) * pii / 180)))
            'Se suman en serie
            K_Washer_Bolt = 1 / (1 / temp1 + 1 / temp2)
        End If
    End If
    '-----------------------------------------
    
    
    '-----------------------------------------
    'Calculo de la rigidez de la arandela de la tuerca
    Dim K_Washer_Nut As Double
    If Nut_Washer_Thickness = 0 Or Nut_Washer_Ext_Dia = 0 Then
        K_Washer_Nut = 10 ^ 100
    Else
        If Nut_Washer_Ext_Dia > WEDN Then
            K_Washer_Nut = K_Tronco_Cono1(Nut_Washer_Young, Alfa, Bore_Dia, Nut_Head_Ext_Dia, Nut_Washer_Thickness)
        ElseIf Nut_Head_Ext_Dia >= WEDN Then
            K_Washer_Nut = K_Cilindro_Hueco(Nut_Washer_Young, Bore_Dia, WEDN, 0, Nut_Washer_Thickness)
        Else
            'Cálculo de la parte cónica
            temp3 = K_Tronco_Cono1(Nut_Washer_Young, Alfa, Bore_Dia, Nut_Head_Ext_Dia, (0.5 * (WEDN - Nut_Head_Ext_Dia)) / (Tan((Alfa) * pii / 180)))
            'Cálculo de la parte cilindrica
            temp4 = K_Cilindro_Hueco(Nut_Washer_Young, Bore_Dia, WEDN, 0, Nut_Washer_Thickness - (0.5 * (WEDN - Nut_Head_Ext_Dia)) / (Tan((Alfa) * pii / 180)))
            'Se suman en serie
            K_Washer_Nut = 1 / (1 / temp3 + 1 / temp4)
        End If
    End If
    '-----------------------------------------
    
    
    '-----------------------------------------
    'Calculo de la rigidez total
    Dim K_Total As Double
    K_Total = 1 / (1 / K_Washer_Bolt + 1 / K_Stack_without_washers + 1 / K_Washer_Nut)
    '-----------------------------------------
    
    
    If selector = "K_Stack_without_washers" Then      'Calculo de la rigidez del stack de piezas sin arandelas
        K_Stack = K_Stack_without_washers
    ElseIf selector = "K_Washer_Bolt" Then                'Calculo de la rigidez de la arandela de la cabeza del fástener
        K_Stack = K_Washer_Bolt
    ElseIf selector = "K_Washer_Nut" Then                 'Calculo de la rigidez de la arandela de la tuerca
        K_Stack = K_Washer_Nut
    ElseIf selector = "K_Total" Then                          'Calculo de la rigidez total
        K_Stack = K_Total
    ElseIf selector = "tabla" Then
        K_Stack = tabla
    End If
    
End Function



