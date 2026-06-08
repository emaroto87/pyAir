Attribute VB_Name = "Ecuaciones"
Function recorrido(fun As String, A As Double, B As Double, error As Double)
    'a es el punto de partida
    'b es el tamaño del paso de resolución
    'error es una variable que no se usa en esta función. Está puesta para que todas las funciones de solver
    ' tengan el mismo número de inputs
    
    Dim error1 As Double
    Dim error2 As Double
    Dim n As Double
    Dim t1 As Double, t2 As Double
            
    t1 = Timer
            
    'se inicializan ambos errores con el mismo valor (el de "a" de partida)
    error1 = eva(fun, CStr(A)) + 1
    error2 = eva(fun, CStr(A))
      
    Do While Abs(error2) < Abs(error1)
        error1 = error2
        A = A + B
        error2 = eva(fun, CStr(A))
        n = n + 1
    Loop
    
    t2 = Timer
    
    Dim vec(3) As Double
    vec(0) = A
    vec(1) = n
    vec(2) = t2 - t1
    recorrido = vec
End Function

Function bolzano(fun As String, A As Double, B As Double, error As Double)
    ' Devuelve el valor de una ecuación que se asume es contínua en el intervalo [a,b]
    
    ' Declaramos las variables necesarias:
    Dim fa As Double
    Dim fb As Double
    Dim c As Double
    Dim fc As Double
    Dim n As Integer
    Dim t1 As Double, t2 As Double
            
    t1 = Timer
            
    ' "a" es el inicio del intervalo donde se desea buscar la solución
    ' "b" es el final del intervalo donde se desea buscar la solución
    ' "error" es la precisión con la que se desea encontrar la solución

    ' En primer lugar verificamos que existe una solución:
    fa = eva(fun, CStr(A))
    fb = eva(fun, CStr(B))
    'fc = eva(fun, CStr(c))      ' Esto era un error ya que "c" se inicializa en cero y fc no tiene porque existir
    fc = error + 1                  ' De esta forma nos aseguramos de que siempre se entra en el bucle
    If fa * fb > 0 Then
        Exit Function
    End If
      
    ' Una vez hecha la comprobación, pasamos a hallar la solución:
    Do While Abs(fc) > error
        c = (A + B) / 2
        fc = eva(fun, CStr(c))
        
        ' Si se encontrara el exacto (poco probable) se da el resultado y salimos
        If fc = 0 Then
            bolzano = c
            Exit Function
        End If
        
        ' Si no se ha encontrado el exacto (muy probable) se calcula en nuevo intervalo
        If fa * fc > 0 Then
            A = c
        Else
            B = c
        End If
                      
        n = n + 1
    Loop
    
    t2 = Timer
    
    Dim vec(3) As Double
    vec(0) = c
    vec(1) = n
    vec(2) = t2 - t1
    bolzano = vec
End Function

Function tangente(fun As String, A As Double, B As Double, error As Double)
    ' Esta método de resolución no siempre converge.
    ' La convergencia depende del punto de partida elegido
    ' Devuelve el valor de una ecuación que se asume es contínua y derivable en la zona de convergencia
    Dim fa As Double
    Dim m As Double
    Dim n As Integer
    Dim t1 As Double, t2 As Double

    t1 = Timer

    ' Se inicializa la variable fa con un valor válido para el algoritmo
    fa = eva(fun, CStr(A))
            
    Do While Abs(fa) > error
        ' Se calcula el valor de la función en el punto de partida
        fa = eva(fun, CStr(A))
        ' Se calcula la pendiente de la recta tangente en la zona de partida
        m = (eva(fun, CStr(A + 0.0000001)) - (fa)) / 0.0000001
        A = A - (fa / m)
        n = n + 1
    Loop
    
    t2 = Timer
    
    Dim vec(3) As Double
    vec(0) = A
    vec(1) = n
    vec(2) = t2 - t1
    tangente = vec
End Function
