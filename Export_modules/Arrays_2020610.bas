Attribute VB_Name = "Arrays_2020610"
Function maxrange(vec As Range)
    ' Returns the max range of an array MAX(rows;cols)
    ' Declaration of variables
    Dim L As Integer
    Dim rows As Integer
    Dim cols As Integer
    
    On Error GoTo check
            
    ' Search for the number or rows and collunms in the array
    rows = UBound(vec(), 1)
    cols = UBound(vec(), 2)
    
    
    ' Get the max between rows and collunms
    If cols > rows Then
        L = cols
    Else
        L = rows
    End If
    
    
' If range is 1x1 UBound returns an error so result is "1" and exit function
check:
    If rows = Empty Then
        maxrange = 1
        Exit Function
    End If
' Else, return the max between rows and cols
    
    
    maxrange = L
End Function

Function sortvec(vec As Range)
    ' Ordena los valores de un rango y los devuelve en un vector
    ' Declaration of variables
    Dim L As Integer
    Dim n As Double
    Dim m As Double
    Dim aux As Double
            
    L = maxrange(vec)
    
    ' Copy the values from the range to the auxiliar vector
    ReDim vec2(L) As Double
    For n = 1 To (L)
        vec2(n - 1) = vec(n)
    Next n
    
    ' Sort vector
    For n = 1 To L
        For m = 0 To (L - 2)
            If vec2(m) > vec2(m + 1) Then
                aux = vec2(m)
                vec2(m) = vec2(m + 1)
                vec2(m + 1) = aux
            End If
        Next m
    Next n

    sortvec = vec2
End Function

Function duplicados(ran As Range)
    ' Elimina los valores duplicados de un rango y devuelve los no duplicados en un vector
    Dim n As Double
    Dim m As Double
    Dim j As Double
    Dim Position As Double
    
    j = 0
    Position = 0
        
    ' Find the length of the range
    Dim L As Integer
    L = maxrange(ran)
    
    ' New vector with le length of the range
    ReDim vec(L, 2) As String
    
    For m = 1 To (L)                    ' Go through the range
        For n = 0 To (L - 1)            ' Go through the vector
            If ran(m) = vec(n, 0) Then
                j = j + 1
            End If
        Next n
        If j = 0 Then
            vec(Position, 0) = ran(m)
            Position = Position + 1
        End If
        j = 0
    Next m

    ' Find how many times each value is repeated
    ReDim vec2(L) As Double
    For n = 0 To (L - 1)                    ' Go through the vector
        For m = 1 To (L)                    ' Go through the range
            If vec(n, 0) = ran(m) Then
                vec2(n) = vec2(n) + 1
            End If
            If vec2(n) <> 0 Then
                vec(n, 1) = CStr(vec2(n))   ' Fill the number of repeats
            End If
        Next m
    Next n

    duplicados = vec
End Function


Function duplicados_vect(entrada() As Double)
    ' Función que elimina los duplicados de un vector y lo devuelve en otro vector
    ' Solo sirve para valores numéricos

    ' Se declaran las variables
    Dim n As Double
    Dim m As Double
    Dim j As Double
    Dim Position As Double
    
    j = 0
    Position = 0
        
    ' Se halla la longitud del vector
    Dim L As Integer
    L = UBound(entrada())

    ' Se declara el nuevo vector
    ReDim salida(L) As Double
    
    For m = 0 To (L)                ' Se recorre el vector de entrada
        For n = 0 To (L)            ' Se recorre el vector de salida
            If entrada(m) = salida(n) Then
                j = j + 1
            End If
        Next n
        If j = 0 Then
            salida(Position) = entrada(m)
            Position = Position + 1
        End If
        j = 0
    Next m

    duplicados_vect = salida
End Function


Function ordena_vect(entrada() As Double)
    ' Ordena los valores de un vector y los devuelve en otro vector
    ' Declaration of variables
    Dim L As Integer
    Dim n As Double
    Dim m As Double
    Dim aux As Double
            
    L = UBound(entrada())
    
    ' Copy the values from the input vector to the output vector
    ReDim salida(L) As Double
        salida = entrada
    
    ' Sort vector
    For n = 0 To L
        For m = 0 To (L - 1)
            If salida(m) > salida(m + 1) Then
                aux = salida(m)
                salida(m) = salida(m + 1)
                salida(m + 1) = aux
            End If
        Next m
    Next n

    ordena_vect = salida
End Function

Sub dependientes()
    ' Muestra las flechas de dependencia de todas las celdas de un rango seleccionado
    Dim m As Double
    Dim n As Double
    Dim p As Double
    n = Selection.Cells.Count
    Application.DisplayStatusBar = True
    
    For Each cell In Selection
        m = m + 1
        p = Round(100 * m / n)
        cell.ShowDependents
        Application.StatusBar = p & "% Cell " & m & " of " & n
    Next
    'Application.DisplayStatusBar = False
End Sub

Sub precedentes()
    ' Muestra las flechas de precedencia de todas las celdas de un rango seleccionado
    Dim m As Double
    Dim n As Double
    Dim p As Double
    n = Selection.Cells.Count
    Application.DisplayStatusBar = True
    
    For Each cell In Selection
        m = m + 1
        p = Round(100 * m / n)
        cell.ShowPrecedents
        Application.StatusBar = p & "% Cell " & m & " of " & n
    Next
    'Application.DisplayStatusBar = False
End Sub
Sub bajar_una_fila()
    ActiveCell.Offset(1, 0).Select
End Sub
Sub modif_columna()
    ActiveCell.Value = ActiveCell.Value & ".els"
    ActiveCell.Offset(1, 0).Select
End Sub
Sub move_rename(carpeta_origen As String, carpeta_destino As String, fichero_origen As String, fichero_destino As String)
    'mover + renombrar archivos de una carpeta a otra
    'carpeta_origen = "D:\Folder"
    'carpeta_destino = "D:\Folder2"
    'fichero_origen = "rrr.dat"
    'fichero_destino = "www.dat"
    
    origen = carpeta_origen & "\" & fichero_origen
    destino = carpeta_destino & "\" & fichero_destino
    
    Name origen As destino

End Sub
Sub copy_paste_rename(carpeta_origen As String, carpeta_destino As String, fichero_origen As String, fichero_destino As String)
    
    origen = carpeta_origen + "\" + fichero_origen
    destino = carpeta_destino + "\" + fichero_destino
    
    Set fss = CreateObject("Scripting.FileSystemObject")
    fss.CopyFile origen, destino

End Sub
Sub Rows_Delete()
    'eliminar filas en blanco(dispuestas una rellena, una en blanco,etc.)
    For i = 3869 To 1 Step -2
        rows(i).Select
        Selection.Delete Shift:=xlUp
    Next
    'Si la haces de arriba hacia abajo al eliminar filas cambia la numeración y se jode
End Sub

Function compactar(ran As Range)
    ' Compacta el contenido de un rango eliminando los huecos vacios
    
    ' Se declaran las variables
    Dim n As Double
    Dim m As Double
    Dim j As Double
    Dim Position As Double
    
    m = 1
        
    ' Find the length of the range
    Dim L As Integer
    L = maxrange(ran)
    
    ' New vector with the length of the range
    ReDim vec(1 To L, 1 To 1) As String
    
    For n = 1 To L
        If ran(n, 1) <> "" Then
            vec(m, 1) = ran(n, 1)
            m = m + 1
        End If
    Next n

    compactar = vec

End Function

Function compactar_horizontal(ran As Range)
    ' Compacta el contenido de un rango eliminando los huecos vacios
    
    ' Se declaran las variables
    Dim n As Double
    Dim m As Double
    Dim j As Double
    Dim Position As Double
    
    m = 1
        
    ' Find the length of the range
    Dim L As Integer
    L = maxrange(ran)
    
    ' New vector with the length of the range
    ReDim vec(1 To 1, 1 To L) As String
    
    For n = 1 To L
        If ran(1, n) <> "" Then
            vec(1, m) = ran(1, n)
            m = m + 1
        End If
    Next n

    compactar_horizontal = vec

End Function


Function BuscaTablaXYZ(X As Double, Y As Double, ranx As Range, rany As Range, ranz As Range)

    Dim n As Double
    n = maxrange(ranx)

    Dim i As Double
    i = 1
    Dim exit_flag As Double
    exit_flag = 0
    
    While exit_flag = 0
        If ranx(i) = X Then
            If rany(i) = Y Then
                exit_flag = 1
            End If
            If i = n Then
                exit_flag = 1
            End If
        End If
        i = i + 1
    Wend
    
    BuscaTablaXYZ = ranz(i - 1)
       

End Function
