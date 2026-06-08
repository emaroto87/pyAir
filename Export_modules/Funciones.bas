Attribute VB_Name = "Funciones"
Function eva(fun As String, X As String)
    Dim S As String
    S = Replace(fun, "x", X)
    eva = Evaluate(S)
End Function

Function integral(fun As String, A As Double, B As Double, steps As Double)
    ' Calcula la integral de una función
    
    ' Declaración de las variables
    Dim delta As Double
    Dim result As Double
    Dim n As Double
    Dim f As Double
        
    ' Cálculo del step de integración
    delta = (B - A) / steps
    result = 0
    
    ' Cálculo de la integral
    For n = 0 To steps
        f = eva(fun, A + n * delta - delta / 2)
        result = result + f * delta
    Next n
    
    integral = result
End Function

Function derivada(fun As String, A As Double, B As Double)
    ' Función que calcula la derivada de una función en un punto
    
    Dim fba As Double
    Dim fab As Double
    Dim m As Double
        
    fba = eva(fun, A - B)
    fab = eva(fun, A + B)
    
    m = (fab - fba) / (2 * B)
    derivada = m
End Function
