Attribute VB_Name = "RemateFinal"
Sub remate()
    Call zoom
    Call comentarios
    Call cella1
    Call borra_nombres
    Worksheets(1).Activate
End Sub

Function numhojas()
    ' Cuenta el número de hojas de calculo que tiene el libro
    ' Si en lugar de "Worksheets.Count" se usa "sheets.Count" devuelve el número total de hojas (no solo las de cálculo)
    numhojas = Worksheets.Count
End Function

Function numcoments()
    ' Cuenta el número de comentarios que hay en la hoja activa
    numcoments = ActiveSheet.Comments.Count
End Function

Sub indice()
    ' Crea un índice de todas las hojas en una columna en una posición fija
    ' Primero Declaramos todas las variables que vamos a necesitar
    Dim nombre As String
    Dim i As Integer
    
    ' Se crea la cabecera del índice
    Cells(1, 10).Activate
    ActiveCell.Value = "ÍNDICE"
    ActiveCell.Font.Bold = True
    ActiveCell.HorizontalAlignment = xlCenter

    ' Recorremos todas las hojas
    For i = 1 To numhojas
        ' Averiguamos el nombre de la hoja
        nombre = Worksheets(i).Name
        ' Activamos la siguiente celda
        Cells(2 + i, 10).Activate
        ' Ponemos el nombre en la celda activa
        ActiveCell.Value = nombre
        ' Cambiamos el color de la celda por el de la "tab" sólamente si ésta tiene color
        If Worksheets(i).Tab.color <> Empty Then
            Selection.Interior.color = Worksheets(i).Tab.color
        Else    ' En el caso de no tener color lo forzamos a blanco ya que para Excel, si no hay color lo pone en negro
            Selection.Interior.color = RGB(255, 255, 255)
        End If
    Next i
End Sub

Sub zoom()
    ' Pone todas las hojas con zoom del 50%
    ' Primero Declaramos todas las variables que vamos a necesitar
    Dim i As Integer
    
    ' Recorremos todas las hojas
    For i = 1 To numhojas
        Worksheets(i).Activate
        ActiveWindow.zoom = 50
    Next i
    Worksheets(1).Activate
End Sub

Sub comentarios()
    Dim i As Integer
    For i = 1 To numhojas
        Worksheets(i).Activate
        Dim cmt As Comment
        For Each cmt In ActiveSheet.Comments
            cmt.Shape.Top = cmt.Parent.Top + 5
            cmt.Shape.Left = cmt.Parent.Offset(0, 1).Left + 5
        Next
    Next i
End Sub

Sub cella1()
    ' Activa la celda A1 en todas las hojas del libro
    Dim i As Integer
    For i = 1 To numhojas
        Worksheets(i).Activate
        Cells(1, 1).Activate
    Next i
    Worksheets(1).Activate
End Sub


Sub objetos()
    For Each hoja In ThisWorkbook.Sheets
        For Each objeto In hoja.DrawingObjects
            objeto.Delete
        Next
    Next
End Sub


Sub borra_nombres()
    ' Subrutina que borra todos los nombres del administrador de nombres que contengan #REF!

    ' Se declaran las variables
    Dim n As Integer
    Dim i As Integer
    Dim j As Integer
    Dim refersto_nombre As String
    Dim checkref As String
    
    j = 1
    ' Se cuenta el número de nombres que hay en el administrador de nombres
    
    Do While j = 1
        j = 0
        n = ActiveWorkbook.Names.Count
        For i = 1 To n
            If j = 0 Then
                refersto_nombre = ActiveWorkbook.Names(i).RefersTo
                checkref = Left(refersto_nombre, 6)
                If checkref = "=#REF!" Then
                    ActiveWorkbook.Names(i).Delete
                    j = 1
                End If
            End If
        Next i
    Loop

End Sub


Sub Links()

    Dim N_Hojas As Double
    N_Hojas = Worksheets.Count
    
    Dim Nombre_Hoja As String
    
        For i = 1 To N_Hojas
            'ruta_archivo = Range("F" & i).value
            'ruta = archivo_con_ruta(ruta_archivo)(1)
            'archivo = archivo_con_ruta(ruta_archivo)(2)
            'If FILE_EXIST(ruta_archivo) Then
                Nombre_Hoja = Range("J" & i + 2).Value
                Range("K" & i + 2).Activate
                ActiveSheet.Hyperlinks.Add Anchor:=Selection, Address:="", SubAddress:= _
                    "'" & Nombre_Hoja & "'!A1", TextToDisplay:="Link"
                    Selection.Font.Size = 11
            'End If
        Next i
End Sub



