Attribute VB_Name = "Vigas"
Function vigas_no_lineal_apoyos(Young As Double, Inercia As Double, longitud As Double, fuerzas As Range, momentos As Range, Fdistribuidas As Range, CContorno As String, pasos As Range, E_interruptor As Range, E_distribucion As Range, I_interruptor As Range, I_distribucion As Range, apoyos As Range, max_iter As Double, Coef_Conv As Double, Axial_Interruptor As Double, Axial_Carga As Double, Axial_Error As Double, Axial_Max_IT As Double, Axial_Interruptor_Flecha_Inicial As Double, Axial_Flecha_Inicial As Range)
    'En principio, es posible que haya que limitar el uso a la condición de contorno "Supported-Supported"
    
    'Se declaran las variables de uso general
    Dim n As Double
    Dim m As Double
    Dim contador As Double
        
    'Se calcula el tamaño de la tabla
    Dim npasos As Double
    npasos = UBound(pasos(), 1)
    
    'Se calcula la viga inicial lineal, con el momento no lineal a cero
    ReDim result_ini(npasos, 6) As Double
    ReDim result_fin(npasos, 6) As Double
    
    'Se crea el vector columna que va a almacenar el histórico de convergencia
    ReDim Historico_Convergencia(npasos, 0) As Double
    
    If Axial_Interruptor <> 1 Then
        Exit Function
    End If
    
    If CContorno <> "Supported-Supported" Then
        Exit Function
    End If
    
    ReDim Momento_No_Lineal_inicial(0 To npasos, 0 To 0) As Double 'Vector columna inicializado en cero para poder pasarlo como argumento
    result_ini = vigaslineal_apoyos(Young, Inercia, longitud, fuerzas, momentos, Fdistribuidas, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, apoyos, max_iter, Coef_Conv, Momento_No_Lineal_inicial)
    
    Dim Max_Flecha_ini As Double
    Dim Max_Flecha_fin As Double
    Dim incremento_flecha As Double
    ReDim Momento_No_Lineal_carga_transversal(npasos, 0) As Double
    Dim Max_Momento_No_Lineal_carga_transversal As Double 'Para check
    incremento_flecha = 1 ' para que entre en el bucle
    
    'Se inicia el bucle
    Do While incremento_flecha > Axial_Error And contador < Axial_Max_IT And contador < npasos 'Se para cuando se alcanza un incremento de flecha inferior al valor introducido por el usuario o bien se alcanzan el número máximo de iteraciones pedido por el usuario o bien se alcanza el número límite de "npasos" iteraciones (aunque el usuario haya pedido más) ya que es el número máximo de interaciones que se puede incluir en la tabla de resultados
        Max_Flecha_ini = 0
        Max_Flecha_fin = 0
        Max_Momento_No_Lineal_carga_transversal = 0 'Para check

        For n = 0 To npasos
            Momento_No_Lineal_carga_transversal(n, 0) = result_ini(n, 4) * Axial_Carga
        Next n
        
        
        '----------------------------------------------------------------------------------------------------------------------------------------------------------------
        'Si se quiere tener en cuenta el desplazamiento lateral inicial forzado:
        If contador = 0 Then
            If Axial_Interruptor_Flecha_Inicial = 1 Then
                For n = 0 To npasos
                    Momento_No_Lineal_carga_transversal(n, 0) = Momento_No_Lineal_carga_transversal(n, 0) + Axial_Flecha_Inicial(n + 1, 1) * Axial_Carga
                Next n
            End If
        End If
        '----------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        For n = 0 To npasos 'Para check
            If Abs(Momento_No_Lineal_carga_transversal(n, 0)) > Max_Momento_No_Lineal_carga_transversal Then
                Max_Momento_No_Lineal_carga_transversal = Abs(Momento_No_Lineal_carga_transversal(n, 0))
            End If
        Next n
        
        
        result_fin = vigaslineal_apoyos(Young, Inercia, longitud, fuerzas, momentos, Fdistribuidas, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, apoyos, max_iter, Coef_Conv, Momento_No_Lineal_carga_transversal)
        
        'Se calcula la máxima flecha inicial y final
        For n = 0 To npasos
            If Abs(result_ini(n, 4)) > Max_Flecha_ini Then
                Max_Flecha_ini = Abs(result_ini(n, 4))
            End If
            If Abs(result_fin(n, 4)) > Max_Flecha_fin Then
                Max_Flecha_fin = Abs(result_fin(n, 4))
            End If
        Next n
        
        incremento_flecha = Max_Flecha_fin / Max_Flecha_ini - 1
        
        'Se copia el resultado inicial al final para la siguiente iteración si ésta fuera necesaria
        For n = 0 To npasos
            For m = 0 To 6
                result_ini(n, m) = result_fin(n, m)
            Next m
        Next n
        
        'Se anota el coeficiente de convergencia
        Historico_Convergencia(contador, 0) = incremento_flecha

        contador = contador + 1
    Loop
    
    'Se copia el vector de convergencia al vector de resultados
    For n = 0 To npasos
        result_fin(n, 6) = Historico_Convergencia(n, 0)
    Next n
    
    vigas_no_lineal_apoyos = result_fin
End Function


Function vigaslineal_apoyos(Young As Double, Inercia As Double, longitud As Double, fuerzas As Range, momentos As Range, Fdistribuidas As Range, CContorno As String, pasos As Range, E_interruptor As Range, E_distribucion As Range, I_interruptor As Range, I_distribucion As Range, apoyos As Range, max_iter As Double, Coef_Conv As Double, Optional Momento_No_Lineal = 0)
    'Esta función es un bucle de control de la función "vigaslineal" que a su vez resuelve vigas de 1 tramo con diferentes condiciones de contorno pero sin apoyos intermedios
    'En esta función realmente no se resuelve la viga.
    'En esta función se busca la convergencia para encontrar las fuerzas que equilibran la viga de forma que actúen como apoyos de la misma
    'Al tratarse de una función "pesada", a lo largo de este código se buscar optimizar para resolver el mínimo número de veces posible la función "vigaslineal"
    
    Dim n As Double
    Dim m As Double
    
    'Se contabilizan en número de fuerzas, momentos, fuerzas distribuidas, pasos y apoyos del problema, estén activos o no
    Dim nfuerzas As Double
    Dim nmomentos As Double
    Dim nfdist As Double
    Dim npasos As Double
    Dim napoyos As Double
    nfuerzas = UBound(fuerzas(), 1)
    nmomentos = UBound(momentos(), 1)
    nfdist = UBound(Fdistribuidas(), 1)
    npasos = UBound(pasos(), 1)
    napoyos = UBound(apoyos(), 1)
    
    'Se pasa el rango de fuerzas a un vector de fuerzas
    'Antes de añadir los apoyos intermedios, las fuerzas se introducían mediante rango como los demás input
    'Para añadir los apoyos intermedios se necesita añadir tantas fuerzas como apoyos existan, por ese motivo hay que pasar la información de un rango a un array que pueda aumentar de tamaño
    ReDim v_fuerzas(1 To nfuerzas, 1 To 3) As Double
    For n = 1 To nfuerzas
        For m = 1 To 3
            v_fuerzas(n, m) = fuerzas(n, m)
        Next m
    Next n
    
    'Se pasan los rangos de momentos y fuerzas distribuidas a respectivos vectores
    'Antes de añadir los apoyos intermedios, los momentos y las fuerzas distribuidas se introducían mediante rango como los demás input
    'Para añadir los apoyos intermedios se necesita calcular la rigidez de la viga en las posiciones de los apoyos. Esto supone cargar la viga individualmente con fuerzas en las posiciones de los apoyos y no puede haber momentos ni fuerzas distribuidas. Por ese motivo, es necesario tener la capacidad de modificar (anular) los momentos o cargas distribuidas que hubiera introducido el usuario (pero solo para estos cálculos de rigidez)
    ReDim v_momentos(1 To nmomentos, 1 To 3) As Double
    For n = 1 To nmomentos
        For m = 1 To 3
            v_momentos(n, m) = momentos(n, m)
        Next m
    Next n
    ReDim v_Fdistribuidas(1 To nfdist, 1 To 5) As Double
    For n = 1 To nfdist
        For m = 1 To 5
            v_Fdistribuidas(n, m) = Fdistribuidas(n, m)
        Next m
    Next n
    
    
    'Se cuentan el número de apoyos introducidos que realmente están activos
    Dim N_apoyos_activos As Double
    N_apoyos_activos = 0
    For n = 1 To nfuerzas
        If apoyos(n, 3) = 1 Then
            N_apoyos_activos = N_apoyos_activos + 1
        End If
    Next n
    
    'Si no hay apoyos se calcula la viga y salimos de la función
    If N_apoyos_activos = 0 Then
        vigaslineal_apoyos = vigaslineal(Young, Inercia, longitud, v_fuerzas, v_momentos, v_Fdistribuidas, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, Momento_No_Lineal)
        Exit Function
    End If
    

    'Si sí hay apoyos pasamos a calcular la rigidez de los mismos
    
    
    '----------------------------------------------------------------------------------------
    'Se calcula la rigidez de la estructura en cada uno de los puntos en los que hay un apoyo
    'Para ello se resuelve la viga para cada apoyo con una carga unitaria en la posición del apoyo
    ReDim rigidez_apoyos(1 To napoyos) As Double
    ReDim posicion_apoyos(1 To napoyos, 1 To 3) As Double
    ReDim fuerzas_calculo_k(1 To 1, 1 To 3) As Double 'Se crea el vector de fuerzas con una única fuerza para el cálculo de la rigidez en el apoyo
    ReDim momentos_calculo_k(1 To 1, 1 To 3) As Double 'Se crea el vector de momentos vacío y con el interrumpor anulado (cero por defecto en variable double)
    ReDim F_dist_calculo_k(1 To 1, 1 To 5) As Double 'Se crea el vector de fuerzas distribuidas vacío y con el interrumpor anulado (cero por defecto en variable double)
    ReDim result_temp(1 To npasos + 1, 6) As Double
    Dim Position_inf As Double
    Dim Position_sup As Double
    Dim Position As Double
    Dim Flecha_calculo_k As Double
    Dim Carga_Unitaria As Double
    Carga_Unitaria = 1
    
    For n = 1 To napoyos
        If apoyos(n, 3) = 1 Then 'Solo se calcula la rigidez si el apoyo está activo
            'Se rellena el vector de fuerzas con una única fuerza para el cálculo de la rigidez en el apoyo
            fuerzas_calculo_k(1, 1) = apoyos(n, 1)   'Posición del apoyo
            fuerzas_calculo_k(1, 2) = Carga_Unitaria 'Carga unitaria
            fuerzas_calculo_k(1, 3) = 1              'Interruptor activo
            'Se resuelve la viga
            result_temp = vigaslineal(Young, Inercia, longitud, fuerzas_calculo_k, momentos_calculo_k, F_dist_calculo_k, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, Momento_No_Lineal)
            'Se busca la posición del apoyo en la tabla
            'Position = Application.WorksheetFunction.Match(apoyos(n, 1), pasos, 0)
            For m = 1 To npasos
                If apoyos(n, 1) = pasos(m, 1) Then
                    Position = m - 1
                End If
                If apoyos(n, 1) < pasos(m, 1) Then
                    Position_sup = m - 1
                    Position_inf = m - 2
                    Exit For
                End If
            Next m
            
            'Se almacena la posición de los apoyos para no perderla y poder usarla más tarde (en el bucle de iteración)
            posicion_apoyos(n, 1) = Position
            posicion_apoyos(n, 2) = Position_sup
            posicion_apoyos(n, 3) = Position_inf
            'Una vez encontrada la posición del apoyo, se busca la flecha que existe en ese punto al haber resuelto
            If Position > 0 Then    'Si la posición del apoyo está en la tabla, se busca su flecha
                Flecha_calculo_k = result_temp(Position, 4)
            End If
            If Position = 0 Then    'Si la posición del apoyo está entre dos valores de la tabla, se interpola su flecha
                Flecha_calculo_k = ((result_temp(Position_sup, 4) - result_temp(Position_inf, 4)) / (pasos(m - 0, 1) - pasos(m - 1, 1))) * (apoyos(n, 1) - pasos(Position_inf, 1)) + result_temp(Position_inf, 4)
            End If
            'Conocida la flecha, se puede calcular la rigidez de la viga en el punto del apoyo
            rigidez_apoyos(n) = Abs(Carga_Unitaria / Flecha_calculo_k) 'En principio debería salir positiva, pero por si acaso se pone el abs()
        End If
        
        'Se reinician las variables para evaluar la rigidez en el siguiente apoyo. Previamente, los valores se han guardado en el vector "posicion_apoyos"
        Position = 0
        Position_sup = 0
        Position_inf = 0
        Flecha_calculo_k = 0
    Next n
    '----------------------------------------------------------------------------------------
    
    
    'Conocida la rigidez de los apoyos se puede pasar al proceso iterativo de equilibrar la viga.
    
    
    '----------------------------------------------------------------------------------------
    'CHECK INICIAL
    'Se resuelve la viga con las cargas reales para un primer check --> Aunque poco probable, podría ocurrir que las condiciones de carga externas impuetas por el usuario cumplan por si solas las condiciones de flecha cero en los apoyos. Es decir, serían unos apoyos con reacción nula porque la viga ya "pasa por ahí" para esas cargas.
    result_temp = vigaslineal(Young, Inercia, longitud, v_fuerzas, v_momentos, v_Fdistribuidas, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, Momento_No_Lineal)
    ReDim flecha_temp(1 To napoyos) As Double
    
    'Se recupera la flecha para los apoyos
    Dim check As Double
    check = 0
    For n = 1 To napoyos
        If apoyos(n, 3) = 1 Then    'Pero solo para los apoyos activos
            'Se recupera el valor de las variables para hacer un copy/paste directo del código de arriba
            Position = posicion_apoyos(n, 1)
            Position_sup = posicion_apoyos(n, 2)
            Position_inf = posicion_apoyos(n, 3)
            'Copy/Paste de arriba:
            If Position > 0 Then    'Si la posición del apoyo está en la tabla, se busca su flecha
                flecha_temp(n) = result_temp(Position, 4)
            End If
            If Position = 0 Then    'Si la posición del apoyo está entre dos valores de la tabla, se interpola su flecha
                flecha_temp(n) = ((result_temp(Position_sup, 4) - result_temp(Position_inf, 4)) / (pasos(m - 0, 1) - pasos(m - 1, 1))) * (apoyos(n, 1) - pasos(Position_inf, 1)) + result_temp(Position_inf, 4)
            End If
            'Fin del Copy/Paste de arriba
            
            'Se chequea si la flecha cumple la condición de apoyo introducida
            If Abs(flecha_temp(n)) < apoyos(n, 2) Then
                check = check + 1
            End If
            
        End If
    Next n
    
    'Si se cumplieran todas las condiciones de sale de la función
    If check = N_apoyos_activos Then
        vigaslineal_apoyos = result_temp
        Exit Function
    End If
    
    'Si no se cumplen, que es lo más probable, se pasa al proceso iterativo
    '----------------------------------------------------------------------------------------
    
    
    
    
    
    '----------------------------------------------------------------------------------------
    'PROCESO ITERATIVO DE BÚSQUEDA DE LA SOLUCIÓN
    Dim N_iter As Double
    'Se crea el vector conjunto de fuerzas + apoyos

    'Se copia el vector de fuerzas al vector conjunto de fuerzas + apoyos
    ReDim v_fuerzas_apoyos(1 To nfuerzas + napoyos, 1 To 3) As Double
    For n = 1 To nfuerzas
        For m = 1 To 3
            v_fuerzas_apoyos(n, m) = v_fuerzas(n, m)
        Next m
    Next n
    
    'Se copia la posición de los apoyos y su estado 1/0 al vector conjunto de fuerzas + apoyos
    For n = 1 To napoyos
        For m = 1 To 3 Step 2
            v_fuerzas_apoyos(n + nfuerzas, m) = apoyos(n, m)
        Next m
    Next n
    
    
    
    
    'COMIENZA EL BUCLE
    ReDim Incremento_Carga(1 To napoyos) As Double
    Do While check < N_apoyos_activos And N_iter < max_iter
        check = 0
        
        'Se resuelve la nueva iteración con el vector conjunto de fuerzas + apoyos
        result_temp = vigaslineal(Young, Inercia, longitud, v_fuerzas_apoyos, v_momentos, v_Fdistribuidas, CContorno, pasos, E_interruptor, E_distribucion, I_interruptor, I_distribucion, Momento_No_Lineal)
        
        'Se recupera la flecha para los apoyos (mismo bucle de más arriba)
        For n = 1 To napoyos
            If apoyos(n, 3) = 1 Then    'Pero solo para los apoyos activos
                'Se recupera el valor de las variables para hacer un copy/paste directo del código de arriba
                Position = posicion_apoyos(n, 1)
                Position_sup = posicion_apoyos(n, 2)
                Position_inf = posicion_apoyos(n, 3)
                If Position > 0 Then    'Si la posición del apoyo está en la tabla, se busca su flecha
                    flecha_temp(n) = result_temp(Position, 4)
                End If
                If Position = 0 Then    'Si la posición del apoyo está entre dos valores de la tabla, se interpola su flecha
                    flecha_temp(n) = ((result_temp(Position_sup, 4) - result_temp(Position_inf, 4)) / (pasos(m - 0, 1) - pasos(m - 1, 1))) * (apoyos(n, 1) - pasos(Position_inf, 1)) + result_temp(Position_inf, 4)
                End If
            End If
        Next n
        
        'Se calcula el Incremente de Carga necesario para equilibrar la viga (F=K·X)
        For n = 1 To napoyos
            Incremento_Carga(n) = rigidez_apoyos(n) * flecha_temp(n) * (-1) / Coef_Conv
        Next n
        
        'Se AÑADE el incremento de carga a la carga existente en el apoyo
        For n = 1 To napoyos
            v_fuerzas_apoyos(n + nfuerzas, 2) = v_fuerzas_apoyos(n + nfuerzas, 2) + Incremento_Carga(n)
        Next n
        
        N_iter = N_iter + 1
    Loop
    '----------------------------------------------------------------------------------------
    


    vigaslineal_apoyos = result_temp
    
End Function

Function vigaslineal(Young As Double, Inercia As Double, longitud As Double, fuerzas() As Double, momentos() As Double, Fdistribuidas() As Double, CContorno As String, pasos As Range, E_interruptor As Range, E_distribucion As Range, I_interruptor As Range, I_distribucion As Range, Optional Momento_No_Lineal = 0)
    
    '30-01-2015 --> Se arregla un bug que impedía poner dos fuerzas en el mismo punto --> solo prevalecía la última de la lista
    '30-01-2015 --> Se arregla un bug que impedía poner dos momentos en el mismo punto --> solo prevalecía el último de la lista
    '30-01-2015 --> Se añade interruptor de eficiencia a las fuerzas
    '30-01-2015 --> Se añade interruptor de eficiencia a los momentos
    '30-01-2015 --> Se añade interruptor de eficiencia a las cargas uniformemente distribuidas

    ' Función que resuelve vigas
    ' Se declaran las variables de uso general
    Dim n As Double
    Dim m As Double
       
    
    ' Se hallan la cantidad de pasos que tiene la viga (en el caso de que se modifique la longitud de la tabla que inicialmente es 1001)
    Dim npasos As Double
    npasos = maxrange(pasos)
    
    
    '----------------------------------------------------------
    'Se rellenan el vector del módulo de Young y la inercia
    ReDim Young_vect(npasos, 0) As Double
    ReDim Inercia_vect(npasos, 0) As Double
    
    For n = 0 To npasos
    
        'Young
        If E_interruptor = "Constante" Then
            Young_vect(n, 0) = Young
        ElseIf E_interruptor = "Variable" Then
            Young_vect(n, 0) = E_distribucion(n + 1, 1)
        End If
        
        'Inercia
        If I_interruptor = "Constante" Then
            Inercia_vect(n, 0) = Inercia
        ElseIf I_interruptor = "Variable" Then
            Inercia_vect(n, 0) = I_distribucion(n + 1, 1)
        End If
        
    Next n
    '----------------------------------------------------------
    
    
    ' Se declara el vector que contendrá todos los resultados en función de la longitud de la tabla. Tiene que tener 6 columnas para contener:
        'Fuerzas aplicadas
        'Momentos aplicados
        'Cortante obtenido
        'Momento obtenido
        'Flecha obtenida
        'Giro obtenido
    ReDim result(npasos, 6) As Double   'En realidad se declara con 7: 0-->6, nos reservamos una por si en el futuro fuera necesaria
    
    
    ' Pasamos las cargas introducidas al vector de cargas (columna 0)
    Dim nfuerzas As Double
    'nfuerzas = maxrange(fuerzas) 'Esto era para cuando las fuerzas se introducían por un rango, antes de añadir los apoyos intermedios
    nfuerzas = UBound(fuerzas(), 1)
    For n = 1 To nfuerzas
        If fuerzas(n, 3) = 1 And fuerzas(n, 1) >= 0 And fuerzas(n, 1) <= longitud Then 'Chequeos de aplicacion: Eficiencia 0/1 y que la fuerza esté en el rango de la viga
            For m = 1 To npasos
                If fuerzas(n, 1) = pasos(m) Then
                    result(m - 1, 0) = result(m - 1, 0) + fuerzas(n, 2)
                ElseIf fuerzas(n, 1) > pasos(m) And fuerzas(n, 1) < pasos(m + 1) Then
                    result(m - 1, 0) = result(m - 1, 0) + (pasos(m + 1) - fuerzas(n, 1)) / (pasos(m + 1) - pasos(m)) * fuerzas(n, 2)
                    result(m - 0, 0) = result(m - 0, 0) + (fuerzas(n, 1) - pasos(m)) / ((pasos(m + 1) - pasos(m))) * fuerzas(n, 2)
                End If
            Next m
        End If
    Next n
    
    
    ' Se pasan los momentos aplicados al vector de momentos (columna 1)
    Dim nmomentos As Double
    'nmomentos = maxrange(momentos) 'Esto era para cuando los momentos se introducían por un rango, antes de añadir los apoyos intermedios
    nmomentos = UBound(momentos(), 1)
    For n = 1 To nmomentos
        If momentos(n, 3) = 1 And momentos(n, 1) >= 0 And momentos(n, 1) <= longitud Then  'Chequeos de aplicacion: Eficiencia 0/1 y que el momento esté en el rango de la viga
            For m = 1 To npasos
                If momentos(n, 1) = pasos(m) Then
                    result(m - 1, 1) = result(m - 1, 1) + momentos(n, 2)
                ElseIf momentos(n, 1) > pasos(m) And momentos(n, 1) < pasos(m + 1) Then
                    result(m - 1, 1) = result(m - 1, 1) + (pasos(m + 1) - momentos(n, 1)) / (pasos(m + 1) - pasos(m)) * momentos(n, 2)
                    result(m - 0, 1) = result(m - 0, 1) + (momentos(n, 1) - pasos(m)) / ((pasos(m + 1) - pasos(m))) * momentos(n, 2)
                End If
            Next m
        End If
    Next n

    
    
    
'    ' ESTE ES EL ANTIGUO ALGORITMO PARA CARGAS RECTANGULARES
'    ' Se pasan las cargas uniformemente distribuidas al vector de cargas (columna 0)
'    Dim nfdist As Double
'    'nfdist = maxrange(Fdistribuidas) 'Esto era para cuando las fuerzas distribuidas se introducían por un rango, antes de añadir los apoyos intermedios
'    nfdist = UBound(Fdistribuidas(), 1)
'    Dim carga As Double
'    Dim PosicionCargaIzq As Double
'    Dim posicionCargaDch As Double
'    Dim CargaPuntual As Double
'
'    For n = 1 To nfdist
'        If Fdistribuidas(n, 4) = 1 Then     ' Eficiencia 0/1
'            carga = (Fdistribuidas(n, 2) - Fdistribuidas(n, 1)) * Fdistribuidas(n, 3)   ' Se halla el total de la carga introducida
'            For m = 1 To npasos
'                If Fdistribuidas(n, 1) >= pasos(m) And Fdistribuidas(n, 1) < pasos(m + 1) Then
'                    PosicionCargaIzq = m
'                End If
'                If Fdistribuidas(n, 2) > pasos(m) And Fdistribuidas(n, 2) <= pasos(m + 1) Then
'                    posicionCargaDch = m + 1
'                End If
'            Next m
'            CargaPuntual = carga / (posicionCargaDch - PosicionCargaIzq)
'            For m = PosicionCargaIzq To posicionCargaDch
'                result(m - 1, 0) = result(m - 1, 0) + CargaPuntual * (npasos) / (npasos + 1)
'            Next m
'        End If
'    Next n
    
    
    
    
    
    ' Se pasan las cargas TRAPEZOIDALES al vector de cargas (columna 0)
    Dim nfdist As Double
    'nfdist = maxrange(Fdistribuidas) 'Esto era para cuando las fuerzas distribuidas se introducían por un rango, antes de añadir los apoyos intermedios
    nfdist = UBound(Fdistribuidas(), 1)
    'Dim carga As Double
    Dim PosicionCargaIzq As Double
    Dim posicionCargaDch As Double
    Dim CargaPuntual As Double
    Dim xini As Double
    Dim xfin As Double
    Dim qini As Double
    Dim qfin As Double
    Dim pendiente As Double
    

    For n = 1 To nfdist
        If Fdistribuidas(n, 5) = 1 And Fdistribuidas(n, 2) > Fdistribuidas(n, 1) And Fdistribuidas(n, 1) >= 0 And Fdistribuidas(n, 2) <= longitud Then ' Chequeos de aplicacion: Eficiencia 0/1, intervalo superiror mayor que el inferior y que toda la carga esté dentro de la viga
            For m = 1 To npasos
                If Fdistribuidas(n, 1) >= pasos(m) And Fdistribuidas(n, 1) < pasos(m + 1) Then
                    PosicionCargaIzq = m
                End If
                If Fdistribuidas(n, 2) > pasos(m) And Fdistribuidas(n, 2) <= pasos(m + 1) Then
                    posicionCargaDch = m + 1
                End If
            Next m
            'Para definir la recta que describe la carga necesitamos la pendiente y al menos un punto (el inicial, por ejemplo)
            xini = Fdistribuidas(n, 1)
            xfin = Fdistribuidas(n, 2)
            qini = Fdistribuidas(n, 3)
            qfin = Fdistribuidas(n, 4)
            pendiente = (qfin - qini) / (xfin - xini)
            
            For m = PosicionCargaIzq To posicionCargaDch - 1
                CargaPuntual = pendiente * (1 / 2) * ((pasos(m + 1) - xini) ^ 2 - (pasos(m) - xini) ^ 2) + qini * (pasos(m + 1) - pasos(m))
                result(m - 1, 0) = result(m - 1, 0) + CargaPuntual
            Next m
        End If
    Next n
    
    
    
    
    
    
    'Se declaran las variables para calcular las integrales para las reacciones en los apoyos
    'Se declaran unas nuevas variables para las reacciones (RA, RB, MA, MB) para, en un principio, no reemplazar el cálculo antiguo de reacciones (Rizq, Rdch, Mizq, Mdch)
    'El cálculo de las antiguas reacciones corresponde a las dsitribuciones constantes de "E" e "I" y, por tanto, no tienen integrales --> Este cálculo consume tan poca cpu que se va a conservar a pesar de que no se va a usar --> Es posible que al final si se acabe usando porque, al consumir menos cpu, es posible que interese usarlo siempre que se pueda, especialmente en cargas distribuidas aplicadas en toda la viga ya que generaría muchas integrales
    Dim RA As Double
    Dim MA As Double
    Dim RB As Double
    Dim MB As Double
    'I1 a I8 corresponden viga biempotrada con carga puntual
    Dim I1 As Double
    Dim I2 As Double
    Dim I3 As Double
    Dim I4 As Double
    Dim I5 As Double
    Dim I6 As Double
    Dim I7 As Double
    Dim I8 As Double
    'I9 a I14 corresponden viga biempotrada con momento
    Dim I9 As Double
    Dim I10 As Double
    Dim I11 As Double
    Dim I12 As Double
    Dim I13 As Double
    Dim I14 As Double
    'I15 a I17 corresponden viga empotrada-apoyada con carga puntual
    Dim I15 As Double
    Dim I16 As Double
    Dim I17 As Double
    'I18 a I20 corresponden viga empotrada-apoyada con momento
    Dim I18 As Double
    Dim I19 As Double
    Dim I20 As Double
    'I21 a I23 corresponden viga empotrada-giro_restringido con carga puntual
    Dim I21 As Double
    Dim I22 As Double
    Dim I23 As Double
    'I24 a I25 corresponden viga empotrada-giro_restringido con momento
    Dim I24 As Double
    Dim I25 As Double

    
    '---------------------------------------------------------------------------------------------
    Select Case CContorno   'En función de las condiciones de contorno se calculan reacciones en los apoyos
        Dim Rizq As Double
        Dim Rdch As Double  'En realidad, las reacciones en el extremo derecho no sirven de nada. Es más para check
        Dim Mizq As Double
        Dim Mdch As Double  'En realidad, las reacciones en el extremo derecho no sirven de nada. Es más para check

        Case "Supported-Supported"  'APOYADO - APOYADO
            Rizq = 0
            Rdch = 0
            Mizq = 0
            Mdch = 0
            For m = 1 To npasos
                Rizq = Rizq + ((longitud - pasos(m)) / longitud) * result(m - 1, 0)
                Rdch = Rdch + (pasos(m) / longitud) * result(m - 1, 0)
            Next m
            Rizq = -Rizq
            Rdch = -Rdch
            Dim sumamomento As Double   ' Se calcula el sumatorio de los momentos aplicados
            sumamomento = 0
            For n = 1 To nmomentos
                If momentos(n, 3) = 1 And momentos(n, 1) >= 0 And momentos(n, 1) <= longitud Then 'Eficiencia 0/1
                    sumamomento = sumamomento + momentos(n, 2)
                End If
            Next n
            Rizq = Rizq + sumamomento / longitud    'Se añade a las reacciones el efecto de los momentos
            Rdch = Rdch - sumamomento / longitud    'Se añade a las reacciones el efecto de los momentos
            
        Case "Fixed-Fixed"          'EMPOTRADO - EMPOTRADO
            Rizq = 0
            Rdch = 0
            Mizq = 0
            Mdch = 0
            For m = 1 To npasos
            
                'Cálculo de las reacciones para "E" e "I" constantes
                If E_interruptor = "Constante" And I_interruptor = "Constante" Then
                    Rizq = Rizq - (result(m - 1, 0) * (longitud ^ 3 - 3 * ((pasos(m)) ^ 2) * longitud + 2 * ((pasos(m)) ^ 3))) / (longitud ^ 3) 'Se contabiliza el efecto de las fuerzas
                    Rizq = Rizq + (6 * result(m - 1, 1) * pasos(m) * (longitud - pasos(m))) / (longitud ^ 3) 'Se contabiliza el efecto de los momentos
                    Mizq = Mizq + (result(m - 1, 0) * pasos(m) * (longitud ^ 2 + (pasos(m)) ^ 2 - 2 * pasos(m) * longitud)) / (longitud ^ 2) 'Se contabiliza el efecto de las fuerzas
                    Mizq = Mizq - (4 * result(m - 1, 1) * pasos(m) * longitud - 3 * result(m - 1, 1) * (pasos(m)) ^ 2 - result(m - 1, 1) * (longitud) ^ 2) / (longitud ^ 2) 'Se contabiliza el efecto de los momentos
                End If
                
                I1 = 0
                I2 = 0
                I3 = 0
                I4 = 0
                I5 = 0
                I6 = 0
                I7 = 0
                I8 = 0
                I9 = 0
                I10 = 0
                I11 = 0
                I12 = 0
                I13 = 0
                I14 = 0
                
                'Cálculo de las reacciones para "E" o "I" variables
                If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
                
                    'FUERZAS
                    If result(m - 1, 0) <> 0 Then   'Solo se calculan las integrales para los puntos de carga (columna 0) no nulos
                        
                        'Cálculo de I1 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I1 = I1 + ((longitud - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I2 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I2 = I2 + result(m - 1, 0) * ((((pasos(m, 1) - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n)))
                        Next n
                        
                        'Cálculo de I3 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I3 = I3 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I4 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I4 = I4 + ((longitud - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I5 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I5 = I5 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I6 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I6 = I6 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I7 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I7 = I7 + result(m - 1, 0) * (((pasos(m, 1) - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n)))
                        Next n
                        
                        'Cálculo de I8 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I8 = I8 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Se calculan las reacciones a partir de las integrales
                        RB = (I7 * (I1 + I4) - I2 * (I6 + I8)) / (((I1 + I4) ^ 2) - (I3 + I5) * (I6 + I8))
                        MB = ((I1 + I4) * RB - I7) / (I6 + I8)
                        RA = RA - (result(m - 1, 0) - RB)
                        MA = MA + result(m - 1, 0) * pasos(m, 1) - longitud * RB + MB
                        
                    End If
                    
                    
                    'MOMENTOS
                    If result(m - 1, 1) <> 0 Then   'Solo se calculan las integrales para los puntos de momento (columna 1) no nulos
                        'Cálculo de I9 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I9 = I9 + ((longitud - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I10 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I10 = I10 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I11 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I11 = I11 + ((longitud - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I12 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I12 = I12 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I13 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I13 = I13 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I14 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I14 = I14 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Se calculan las reacciones a partir de las integrales
                        RB = result(m - 1, 1) * ((I13 * (I9 + I11) - I9 * (I13 + I14)) / (-((I9 + I11) ^ 2) + (I10 + I12) * (I13 + I14)))
                        MB = (-(I9 + I11) * RB - I13 * result(m - 1, 1)) / (I13 + I14)
                        RA = RA - RB
                        MA = MA + longitud * RB + MB + result(m - 1, 1)
    
                    End If
                    
                    'Se cambia el nombre a las variables
                    Rizq = RA
                    Mizq = MA
                End If
                
            Next m
            
        Case "Fixed-Supported"      'EMPOTRADO - APOYADO
            Rizq = 0
            Rdch = 0
            Mizq = 0
            Mdch = 0
            For m = 1 To npasos
                'Cálculo de las reacciones para "E" e "I" constantes
                If E_interruptor = "Constante" And I_interruptor = "Constante" Then
                    Rizq = Rizq - (result(m - 1, 0) * (2 * longitud ^ 3 - 3 * longitud * (pasos(m)) ^ 2 + pasos(m) ^ 3)) / (2 * longitud ^ 3) 'Se contabiliza el efecto de las fuerzas
                    Rizq = Rizq + (3 * result(m - 1, 1) * pasos(m) * (2 * longitud - pasos(m))) / (2 * longitud ^ 3) 'Se contabiliza el efecto de los momentos
                    Mizq = Mizq + (result(m - 1, 0) * pasos(m) * (2 * longitud ^ 2 - 3 * longitud * pasos(m) + (pasos(m)) ^ 2)) / (2 * longitud ^ 2) 'Se contabiliza el efecto de las fuerzas
                    Mizq = Mizq - (result(m - 1, 1) * (6 * pasos(m) * longitud - 3 * (pasos(m)) ^ 2 - 2 * longitud ^ 2)) / (2 * longitud ^ 2) 'Se contabiliza el efecto de los momentos
                End If
                
                I15 = 0
                I16 = 0
                I17 = 0
                I18 = 0
                I19 = 0
                I20 = 0
                
                'Cálculo de las reacciones para "E" o "I" variables
                If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
                    
                    'FUERZAS
                    If result(m - 1, 0) <> 0 Then   'Solo se calculan las integrales para los puntos de carga (columna 0) no nulos
                    
                        'Cálculo de I15 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I15 = I15 + (((pasos(m, 1) - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I16 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I16 = I16 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I17 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I17 = I17 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Se calculan las reacciones a partir de las integrales
                        RB = (result(m - 1, 0) * I15) / (I16 + I17)
                        RA = RA - (result(m - 1, 0) - RB) 'Se cambia el signo para ajustar con el criterio de signos que tenía el programa hasta esta versión
                        MA = MA + result(m - 1, 0) * pasos(m, 1) - longitud * RB
                        
                    End If
                    
                    
                    'MOMENTOS
                    If result(m - 1, 1) <> 0 Then   'Solo se calculan las integrales para los puntos de momento (columna 1) no nulos
                        
                        'Cálculo de I18 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I18 = I18 + ((longitud - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I19 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I19 = I19 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I20 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I20 = I20 + (((longitud - pasos(n, 1)) * (longitud - pasos(n, 1))) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Se calculan las reacciones a partir de las integrales
                        RB = (-result(m - 1, 1) * I18) / (I19 + I20)
                        RA = RA - RB
                        MA = MA - (-result(m - 1, 1) - longitud * RB) 'Se cambia el signo para ajustar con el criterio de signos que tenía el programa hasta esta versión
                    End If
                    
                    'Se cambia el nombre a las variables
                    Rizq = RA
                    Mizq = MA
                    
                End If
                
            Next m
        Case "Fixed-Rotation_Fixed" 'EMPOTRADO - GIRO RESTRINGIDO (pero no el desplazamiento vertical)
            Rizq = 0
            Rdch = 0
            Mizq = 0
            Mdch = 0
            For m = 1 To npasos
                'Cálculo de las reacciones para "E" e "I" constantes
                If E_interruptor = "Constante" And I_interruptor = "Constante" Then
                    Rizq = Rizq - result(m - 1, 0) 'Se contabiliza el efecto de las fuerzas
                    Rizq = Rizq + 0 'Se contabiliza el efecto de los momentos
                    Mizq = Mizq + (result(m - 1, 0) * pasos(m) * (2 * longitud - pasos(m))) / (2 * longitud) 'Se contabiliza el efecto de las fuerzas
                    Mizq = Mizq - (result(m - 1, 1) * (pasos(m) - longitud)) / (longitud) 'Se contabiliza el efecto de los momentos
                End If
                
                I21 = 0
                I22 = 0
                I23 = 0
                I24 = 0
                I25 = 0
                
                'Cálculo de las reacciones para "E" o "I" variables
                If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
                    
                    'FUERZAS
                    If result(m - 1, 0) <> 0 Then   'Solo se calculan las integrales para los puntos de carga (columna 0) no nulos
                    
                        'Cálculo de I21 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I21 = I21 + ((pasos(m, 1) - pasos(n, 1)) / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I22 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I22 = I22 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I23 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I23 = I23 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Se calculan las reacciones a partir de las integrales
                        MB = (-result(m - 1, 0) * I21) / (I22 + I23)
                        RA = RA - result(m - 1, 0)
                        MA = MA + result(m - 1, 0) * pasos(m, 1) + MB
                        
                    End If
                    
                    'MOMENTOS
                    If result(m - 1, 1) <> 0 Then   'Solo se calculan las integrales para los puntos de momento (columna 1) no nulos
                    
                        'Cálculo de I24 (la integral va desde el inicio hasta el punto de aplicación de la carga)
                        For n = 1 To m
                            I24 = I24 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                        
                        'Cálculo de I25 (la integral va desde el punto de aplicación de la carga hasta el final)
                        For n = m To npasos - 1
                            I25 = I25 + (1 / (Young_vect(n - 1, 0) * Inercia_vect(n - 1, 0))) * (pasos(n + 1) - pasos(n))
                        Next n
                    
                        'Se calculan las reacciones a partir de las integrales
                        MB = (result(m - 1, 1) * I24) / (I24 + I25)
                        RA = RA + 0
                        MA = MA - (MB - result(m - 1, 1))
                    
                    End If
                        
                    'Se cambia el nombre a las variables
                    Rizq = RA
                    Mizq = MA

                End If

            Next m
            
        Case "Fixed-Free"           'EMPOTRADO - LIBRE
            Rizq = 0
            Rdch = 0
            Mizq = 0
            Mdch = 0
            For m = 1 To npasos
                Rizq = Rizq - result(m - 1, 0)              'Para esta configuración los momentos aplicados no generan cortante
                Mizq = Mizq + result(m - 1, 0) * pasos(m)   'Se contabiliza el efecto de las cargas aplicadas
                Mizq = Mizq + result(m - 1, 1)              'Se contabiliza el efecto de los momentos aplicados
            Next m

    End Select
    
    
    '---------------------------------------------------------------------------------------------
    ' Se rellena el vector del diagrama de esfuerzos cortantes (columna 2)
    result(0, 2) = Rizq + result(0, 0)
    For m = 2 To (npasos)
        result(m - 1, 2) = result(m - 2, 2) + result(m - 1, 0)
    Next m
        'result(npasos - 1, 2) = 0   'La ultima celda del vector la dejamos vacía
    '---------------------------------------------------------------------------------------------
        
        
    '---------------------------------------------------------------------------------------------
    ' Se rellena el vector del diagrama de momentos flectores (columna 3)
    result(0, 3) = Mizq - result(0, 1) 'No hay momento en un extremo apoyado a menos que se aplique directamente
    For m = 2 To (npasos)
        result(m - 1, 3) = result(m - 2, 3) + (pasos(m) - pasos(m - 1)) * (result(m - 2, 2)) - result(m - 1, 1)
    Next m
    ' Se cambia el signo al diagrama de momentos flectores para que tienda a parecerse más a la deformada
    For m = 1 To (npasos)
        result(m - 1, 3) = -result(m - 1, 3)
    Next m
    '---------------------------------------------------------------------------------------------
    
    
    '---------------------------------------------------------------------------------------------
    'Se añade el momento no lineal al momento flector que se acaba de calcular
    On Error GoTo err
    If Momento_No_Lineal = 0 Then   'Si no se ha introducido el momento no lineal como argumento, se crea un argumento inicializado en cero para poder operar con él
        ReDim Momento_No_Lineal(0 To npasos, 0) As Double
    End If
err:
    For m = 1 To (npasos)
        result(m - 1, 3) = result(m - 1, 3) - Momento_No_Lineal(m - 1, 0)
    Next m
    '---------------------------------------------------------------------------------------------

    
    
    Select Case CContorno   'En función de las condiciones de contorno se calcula el giro en el extremo izquierdo
        Dim integral As Double
        Dim giroizq As Double
        Case "Supported-Supported"  'APOYADO - APOYADO
            integral = 0
            For m = 1 To (npasos - 1)
                If E_interruptor = "Constante" And I_interruptor = "Constante" Then
                    integral = integral + ((result(m, 3) + result(m - 1, 3)) / 2) * (pasos(m + 1) - pasos(m)) * (longitud - pasos(m + 1))
                End If
                If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
                    integral = integral + (((result(m, 3) + result(m - 1, 3)) / 2) * (pasos(m + 1) - pasos(m)) * (longitud - pasos(m + 1))) / (Young_vect(m - 1, 0) * Inercia_vect(m - 1, 0))
                End If
            Next m
            If E_interruptor = "Constante" And I_interruptor = "Constante" Then
                giroizq = (1 / (Young * Inercia * longitud)) * integral
            End If
            If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
                giroizq = (1 / longitud) * integral
            End If
        Case "Fixed-Fixed"          'EMPOTRADO - EMPOTRADO
            giroizq = 0
        Case "Fixed-Supported"      'EMPOTRADO - APOYADO
            giroizq = 0
        Case "Fixed-Rotation_Fixed" 'EMPOTRADO - GIRO RESTRINGIDO
            giroizq = 0
        Case "Fixed-Free"           'EMPOTRADO - LIBRE
            giroizq = 0
    End Select
   
   
    '---------------------------------------------------------------------------------------------
    'Se rellena el vector del giro (columna 5)
    result(0, 5) = giroizq
    For n = 2 To npasos
        If E_interruptor = "Constante" And I_interruptor = "Constante" Then
            result(n - 1, 5) = result(n - 2, 5) - (1 / (Young * Inercia)) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1))
        End If
        If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
            result(n - 1, 5) = result(n - 2, 5) - (1 / (Young_vect(n - 2, 0) * Inercia_vect(n - 2, 0))) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1))
        End If
    Next n
    '---------------------------------------------------------------------------------------------



    '---------------------------------------------------------------------------------------------
    ' Se rellena el vector de la flecha (columna 4)
    result(0, 4) = 0    ' En el extremo izquierdo la flecha es cero (para apoyado y empotrado)
    For n = 2 To npasos
        If E_interruptor = "Constante" And I_interruptor = "Constante" Then
            result(n - 1, 4) = result(n - 2, 4) + result(n - 2, 5) * (pasos(n) - pasos(n - 1)) + (1 / (Young * Inercia)) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1)) * (1 / 2) * (pasos(n) - pasos(n - 1))
        End If
        If E_interruptor = "Variable" Or I_interruptor = "Variable" Then
            result(n - 1, 4) = result(n - 2, 4) + result(n - 2, 5) * (pasos(n) - pasos(n - 1)) + (1 / (Young_vect(n - 2, 0) * Inercia_vect(n - 2, 0))) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1)) * (1 / 2) * (pasos(n) - pasos(n - 1))
        End If
    Next n
    '---------------------------------------------------------------------------------------------
    
    vigaslineal = result
    
End Function



'Function vigaslinealmodificada(Young As Double, Inercia As Double, longitud As Double, fuerzas As Range, momentos As Range, Fdistribuidas As Range, CContorno As String, pasos As Range, momentosnolineal() As Double)
'    'PENDIENTE IMPLEMENTAR LAS EFICIENCIAS 0/1
'    ' Es la misma función para resolver vigas lineales pero que acepta un vector con los momentos secundarios debido a la carga axial
'    ' Se declaran las variables de uso general
'    Dim n As Double
'    Dim m As Double
'
'
'    ' Se hallan la cantidad de pasos que tiene la viga (en el caso de que se modifique la longitud de la tabla que inicialmente es 1001)
'    Dim npasos As Double
'    npasos = maxrange(pasos)
'
'
'    ' Se declara el vector que contendrá todos los resultados en función de la longitud de la tabla. Tiene que tener 6 columnas para contener:
'        'Fuerzas aplicadas
'        'Momentos aplicados
'        'Cortante obtenido
'        'Momento obtenido
'        'Flecha obtenida
'        'Giro obtenido
'    ReDim result(npasos, 6) As Double
'
'
'    ' Pasamos las cargas introducidas al vector de cargas (columna 0)
'    Dim nfuerzas As Double
'    nfuerzas = maxrange(fuerzas)
'    For n = 1 To nfuerzas
'        For m = 1 To npasos
'            If fuerzas(n, 1) = pasos(m) Then
'                result(m - 1, 0) = fuerzas(n, 2)
'            ElseIf fuerzas(n, 1) > pasos(m) And fuerzas(n, 1) < pasos(m + 1) Then
'                result(m - 1, 0) = (pasos(m + 1) - fuerzas(n, 1)) / (pasos(m + 1) - pasos(m)) * fuerzas(n, 2)
'                result(m - 0, 0) = (fuerzas(n, 1) - pasos(m)) / ((pasos(m + 1) - pasos(m))) * fuerzas(n, 2)
'            End If
'        Next m
'    Next n
'
'
'    ' Se pasan los momentos aplicados al vector de momentos (columna 1)
'    Dim nmomentos As Double
'    nmomentos = maxrange(momentos)
'    For n = 1 To nmomentos
'        For m = 1 To npasos
'            If momentos(n, 1) = pasos(m) Then
'                result(m - 1, 1) = momentos(n, 2)
'            ElseIf momentos(n, 1) > pasos(m) And momentos(n, 1) < pasos(m + 1) Then
'                result(m - 1, 1) = (pasos(m + 1) - momentos(n, 1)) / (pasos(m + 1) - pasos(m)) * momentos(n, 2)
'                result(m - 0, 1) = (momentos(n, 1) - pasos(m)) / ((pasos(m + 1) - pasos(m))) * momentos(n, 2)
'            End If
'        Next m
'    Next n
'
'
'
'    ' Se pasan las cargas uniformemente distribuidas al vector de cargas (columna 0)
'    Dim nfdist As Double
'    nfdist = maxrange(Fdistribuidas)
'    Dim carga As Double
'    Dim PosicionCargaIzq As Double
'    Dim posicionCargaDch As Double
'    Dim CargaPuntual As Double
'
'    For n = 1 To nfdist
'        carga = (Fdistribuidas(n, 2) - Fdistribuidas(n, 1)) * Fdistribuidas(n, 3)   ' Se halla el total de la carga introducida
'        For m = 1 To npasos
'            If Fdistribuidas(n, 1) >= pasos(m) And Fdistribuidas(n, 1) < pasos(m + 1) Then
'                PosicionCargaIzq = m
'            End If
'            If Fdistribuidas(n, 2) > pasos(m) And Fdistribuidas(n, 2) <= pasos(m + 1) Then
'                posicionCargaDch = m + 1
'            End If
'        Next m
'        CargaPuntual = carga / (posicionCargaDch - PosicionCargaIzq)
'        For m = PosicionCargaIzq To posicionCargaDch
'            result(m - 1, 0) = result(m - 1, 0) + CargaPuntual * (npasos) / (npasos + 1)
'        Next m
'    Next n
'
'    '---------------------------------------------------------------------------------------------
'    Select Case CContorno   'En función de las condiciones de contorno se calculan reacciones en los apoyos
'        Dim Rizq As Double
'        Dim Rdch As Double  'En realidad, las reacciones en el extremo derecho no sirven de nada. Es más para check
'        Dim Mizq As Double
'        Dim Mdch As Double  'En realidad, las reacciones en el extremo derecho no sirven de nada. Es más para check
'
'        Case "Supported-Supported"  'APOYADO - APOYADO
'            Rizq = 0
'            Rdch = 0
'            Mizq = 0
'            Mdch = 0
'            For m = 1 To npasos
'                Rizq = Rizq + ((longitud - pasos(m)) / longitud) * result(m - 1, 0)
'                Rdch = Rdch + (pasos(m) / longitud) * result(m - 1, 0)
'            Next m
'            Rizq = -Rizq
'            Rdch = -Rdch
'            Dim sumamomento As Double   ' Se calcula el sumatorio de los momentos aplicados
'            sumamomento = 0
'            For n = 1 To nmomentos
'                sumamomento = sumamomento + momentos(n, 2)
'            Next n
'            Rizq = Rizq + sumamomento / longitud    'Se añade a las reacciones el efecto de los momentos
'            Rdch = Rdch - sumamomento / longitud    'Se añade a las reacciones el efecto de los momentos
'        Case "Fixed-Fixed"          'EMPOTRADO - EMPOTRADO
'            Rizq = 0
'            Rdch = 0
'            Mizq = 0
'            Mdch = 0
'            For m = 1 To npasos
'                Rizq = Rizq - (result(m - 1, 0) * (longitud ^ 3 - 3 * ((pasos(m)) ^ 2) * longitud + 2 * ((pasos(m)) ^ 3))) / (longitud ^ 3) 'Se contabiliza el efecto de las fuerzas
'                Rizq = Rizq + (6 * result(m - 1, 1) * pasos(m) * (longitud - pasos(m))) / (longitud ^ 3) 'Se contabiliza el efecto de los momentos
'                Mizq = Mizq + (result(m - 1, 0) * pasos(m) * (longitud ^ 2 + (pasos(m)) ^ 2 - 2 * pasos(m) * longitud)) / (longitud ^ 2) 'Se contabiliza el efecto de las fuerzas
'                Mizq = Mizq - (4 * result(m - 1, 1) * pasos(m) * longitud - 3 * result(m - 1, 1) * (pasos(m)) ^ 2 - result(m - 1, 1) * (longitud) ^ 2) / (longitud ^ 2) 'Se contabiliza el efecto de los momentos
'            Next m
'        Case "Fixed-Supported"      'EMPOTRADO - APOYADO
'            Rizq = 0
'            Rdch = 0
'            Mizq = 0
'            Mdch = 0
'            For m = 1 To npasos
'                Rizq = Rizq - (result(m - 1, 0) * (2 * longitud ^ 3 - 3 * longitud * (pasos(m)) ^ 2 + pasos(m) ^ 3)) / (2 * longitud ^ 3) 'Se contabiliza el efecto de las fuerzas
'                Rizq = Rizq + (3 * result(m - 1, 1) * pasos(m) * (2 * longitud - pasos(m))) / (2 * longitud ^ 3) 'Se contabiliza el efecto de los momentos
'                Mizq = Mizq + (result(m - 1, 0) * pasos(m) * (2 * longitud ^ 2 - 3 * longitud * pasos(m) + (pasos(m)) ^ 2)) / (2 * longitud ^ 2) 'Se contabiliza el efecto de las fuerzas
'                Mizq = Mizq - (result(m - 1, 1) * (6 * pasos(m) * longitud - 3 * (pasos(m)) ^ 2 - 2 * longitud ^ 2)) / (2 * longitud ^ 2) 'Se contabiliza el efecto de los momentos
'            Next m
'        Case "Fixed-Free"           'EMPOTRADO - LIBRE
'            Rizq = 0
'            Rdch = 0
'            Mizq = 0
'            Mdch = 0
'            For m = 1 To npasos
'                Rizq = Rizq - result(m - 1, 0)              'Para esta configuración los momentos aplicados no generan cortante
'                Mizq = Mizq + result(m - 1, 0) * pasos(m)   'Se contabiliza el efecto de las cargas aplicadas
'                Mizq = Mizq + result(m - 1, 1)              'Se contabiliza el efecto de los momentos aplicados
'            Next m
'
'    End Select
'
'
'    '---------------------------------------------------------------------------------------------
'    ' Se rellena el vector del diagrama de esfuerzos cortantes (columna 2)
'    result(0, 2) = Rizq + result(0, 0)
'    For m = 2 To (npasos)
'        result(m - 1, 2) = result(m - 2, 2) + result(m - 1, 0)
'    Next m
'        'result(npasos - 1, 2) = 0   'La ultima celda del vector la dejamos vacía
'    '---------------------------------------------------------------------------------------------
'
'
'    ' Se rellena el vector del diagrama de momentos flectores (columna 3)
'    result(0, 3) = Mizq - result(0, 1) 'No hay momento en un extremo apoyado a menos que se aplique directamente
'    For m = 2 To (npasos)
'        result(m - 1, 3) = result(m - 2, 3) + (pasos(m) - pasos(m - 1)) * (result(m - 2, 2)) - result(m - 1, 1)
'    Next m
'    ' Se cambia el signo al diagrama de momentos flectores para que tienda a parecerse más a la deformada
'    For m = 1 To (npasos)
'        result(m - 1, 3) = -result(m - 1, 3)
'    Next m
'        'result(npasos - 1, 3) = 0   'La ultima celda del vector la dejamos vacía
'
'
'
'    ' NOVEDAD CON RESPECTO A LA FUNCIÓN ESTANDARD
'    '___________________________________________________________________
'    'Se añaden al vector de la distribución de momentos (columna 3) los momentos debidos a la carga axial
'    For n = 1 To npasos
'        result(n - 1, 3) = result(n - 1, 3) + momentosnolineal(n - 1, 0)
'    Next n
'    '___________________________________________________________________
'    '---------------------------------------------------------------------------------------------
'
'    Select Case CContorno   'En función de las condiciones de contorno se calcula el giro en el extremo izquierdo
'        Dim integral As Double
'        Dim giroizq As Double
'        Case "Supported-Supported"  'APOYADO - APOYADO
'            integral = 0
'            For m = 1 To (npasos - 1)
'                integral = integral + ((result(m, 3) + result(m - 1, 3)) / 2) * (pasos(m + 1) - pasos(m)) * (longitud - pasos(m + 1))
'            Next m
'            giroizq = (1 / (Young * Inercia * longitud)) * integral
'        Case "Fixed-Fixed"          'EMPOTRADO - EMPOTRADO
'            giroizq = 0
'        Case "Fixed-Supported"      'EMPOTRADO - APOYADO
'            giroizq = 0
'        Case "Fixed-Free"           'EMPOTRADO - LIBRE
'            giroizq = 0
'    End Select
'
'    '---------------------------------------------------------------------------------------------
'    'Se rellena el vector del giro (columna 5)
'    result(0, 5) = giroizq
'    For n = 2 To npasos
'        result(n - 1, 5) = result(n - 2, 5) - (1 / (Young * Inercia)) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1))
'    Next n
'
'    '---------------------------------------------------------------------------------------------
'
'
'    ' Se rellena el vector de la flecha (columna 4)
'    result(0, 4) = 0    ' En el extremo izquierdo la flecha es cero (para apoyado y empotrado)
'    For n = 2 To npasos
'        result(n - 1, 4) = result(n - 2, 4) + result(n - 2, 5) * (pasos(n) - pasos(n - 1)) + (1 / (Young * Inercia)) * ((result(n - 1, 3) + result(n - 2, 3)) / 2) * (pasos(n) - pasos(n - 1)) * (1 / 2) * (pasos(n) - pasos(n - 1))
'    Next n
'
'    vigaslinealmodificada = result
'
'End Function
'
'
