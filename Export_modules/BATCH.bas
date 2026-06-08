Attribute VB_Name = "BATCH"
Sub batch()
    ' Subrutina que loopea casos de carga y construye el envelope (opcional)

    ' Se declaran las variables que se van a usar:
    Dim lc_sheet As String
    Dim lc_batch As String
    Dim lc_label_col As String
    Dim lc_num As Double
    Dim lc_result As String
    Dim refresh As String
    Dim Conservar As String
    Dim envelope_switch As String
    Dim envelope_nombre As String
    Dim Conservar_envelope As String
    Dim Auto_Save As String
    
    ' Se rellenan las variables que se van a usar
    Range("E15").Activate
    lc_sheet = ActiveCell.Value
    
    Range("E16").Activate
    lc_batch = ActiveCell.Value
    
    Range("E17").Activate
    lc_label_col = ActiveCell.Value
    
    Range("E18").Activate
    lc_num = ActiveCell.Value
    
    Range("E19").Activate
    lc_result = ActiveCell.Value

    Range("E20").Activate
    refresh = ActiveCell.Value

    Range("E21").Activate
    Conservar = ActiveCell.Value

    Range("E22").Activate
    envelope_switch = ActiveCell.Value

    Range("E23").Activate
    envelope_nombre = ActiveCell.Value

    Range("E24").Activate
    Conservar_envelope = ActiveCell.Value
    
    Range("E25").Activate
    Auto_Save = ActiveCell.Value


    ' Se desactiva el refresto de pantalla (si es requerido)
    If refresh = "NO" Then
        Application.ScreenUpdating = False
    End If


    ' Se cuenta el número total de hojas que hay en el libro
    Dim Num_Hojas As Integer
    Dim m As Integer
    Dim Check_Hoja_Resultados As Integer
    Num_Hojas = Worksheets.Count
    ' Se comprueba si existe la hoja de resultados
    For m = 1 To Num_Hojas
        If lc_result = Worksheets(m).Name Then
            Check_Hoja_Resultados = 1
        End If
    Next m
    ' Si la hoja no existe se crea
    If Check_Hoja_Resultados <> 1 Then
        Worksheets.Add.Name = lc_result
    End If
     
    ' Se guardan los resultados anteriores
    If Check_Hoja_Resultados = 1 And Conservar = "SI" Then
        Worksheets(lc_result).Copy After:=ActiveWorkbook.Sheets(lc_result)
    End If
     
 
    ' Se copia la columna de las labels de los casos de carga a la hoja de loop
    Worksheets(lc_sheet).Activate
    Range(lc_label_col & ":" & lc_label_col).Copy
    Worksheets(lc_batch).Activate
    Range("A:A").PasteSpecial xlPasteAll
    Application.CutCopyMode = False
    ' Se copia la columna de las labels de los resultados a la hoja de resultados
    Range("J:J").Copy
    Worksheets(lc_result).Activate
    Range("A:A").PasteSpecial xlPasteAll
    Application.CutCopyMode = False
    
' INICIO DE LA LOOP
' -----------------
    Dim n As Double
    Dim Num_Col_Label As Double
    Worksheets(lc_sheet).Activate
    Range(lc_label_col & "1").Activate
    Num_Col_Label = ActiveCell.column
    For n = 1 To lc_num
        Application.DisplayStatusBar = True
        Application.StatusBar = "Load Case " & n & " of " & lc_num
        Worksheets(lc_sheet).Activate
        Range(Columns(Num_Col_Label + n), Columns(Num_Col_Label + n)).Copy
        Worksheets(lc_batch).Activate
        Range("B:B").PasteSpecial xlPasteAll
        Application.CutCopyMode = False
        Application.CalculateFull
        Range("K:K").Copy
        Worksheets(lc_result).Activate
        Range(Columns(n + 1), Columns(n + 1)).PasteSpecial xlPasteAll
        Range(Columns(n + 1), Columns(n + 1)).PasteSpecial xlPasteValues
        Application.CutCopyMode = False
    Next n
' -----------------
' FIN DE LA LOOP
    Application.DisplayStatusBar = False
    If refresh = "NO" Then
        Application.ScreenUpdating = False
    End If

' CONSTRUCCION DEL ENVELOPE (SI ES REQUERIDO)
' -----------------
    If envelope_switch = "SI" Then
        ' Se comprueba si existe la hoja de envelope
        Dim Check_Hoja_envelope As Integer
        For m = 1 To Num_Hojas
            If envelope_nombre = Worksheets(m).Name Then
                Check_Hoja_envelope = 1
            End If
        Next m
        ' Si la hoja no existe se crea
        If Check_Hoja_envelope <> 1 Then
            Worksheets.Add.Name = envelope_nombre
        End If

        ' Se guardan los envelope anteriores
        If Check_Hoja_envelope = 1 And Conservar_envelope = "SI" Then
            Worksheets(envelope_nombre).Copy After:=ActiveWorkbook.Sheets(envelope_nombre)
        End If

        ' Se copia la columna de las labels de los resultados a la hoja de envelope
        Worksheets(lc_batch).Activate
        Range("J:J").Copy
        Worksheets(envelope_nombre).Activate
        Range("A:A").PasteSpecial xlPasteAll
        Application.CutCopyMode = False
        
        ' Se construye la cabecera del envelope
        Range("B1").Activate
        ActiveCell.Value = "MAX"
        Range("C1").Activate
        ActiveCell.Value = "MAX LOAD CASE NUMBER"
        Range("D1").Activate
        ActiveCell.Value = "MAX LOAD CASE"
        Range("E1").Activate
        ActiveCell.Value = "MIN"
        Range("F1").Activate
        ActiveCell.Value = "MIN LOAD CASE NUMBER"
        Range("G1").Activate
        ActiveCell.Value = "MIN LOAD CASE"

        
        ' Se recorren los resultados y se crea el envelope
        Worksheets(lc_batch).Activate
        Range("L2").Activate
        n = 2
        Dim maximo As Double
        Dim max_lc_posicion As Integer
        Dim max_lc As String
        Dim minimo As Double
        Dim min_lc_posicion As Integer
        Dim min_lc As String
        Do While ActiveCell.Value <> Empty
            If ActiveCell.Value = "SI" Then
                Worksheets(lc_result).Activate
                ' Se calcula el máximo del envelope
                maximo = Application.WorksheetFunction.MAX(Range(n & ":" & n))
                max_lc_posicion = Application.WorksheetFunction.Match(maximo, Range(n & ":" & n), 0)
                max_lc = Application.WorksheetFunction.Index(Range("1:1"), max_lc_posicion)
                
                ' Se calcula el mínmo del envelope
                minimo = Application.WorksheetFunction.MIN(Range(n & ":" & n))
                min_lc_posicion = Application.WorksheetFunction.Match(minimo, Range(n & ":" & n), 0)
                min_lc = Application.WorksheetFunction.Index(Range("1:1"), min_lc_posicion)

                Worksheets(envelope_nombre).Activate
                ' Se graban los datos del máximo
                Range("B" & n).Activate
                ActiveCell.Value = maximo
                Range("C" & n).Activate
                ActiveCell.Value = max_lc_posicion - 1
                Range("D" & n).Activate
                ActiveCell.Value = max_lc
                
                ' Se graban los datos del mínimo
                Range("E" & n).Activate
                ActiveCell.Value = minimo
                Range("F" & n).Activate
                ActiveCell.Value = min_lc_posicion - 1
                Range("G" & n).Activate
                ActiveCell.Value = min_lc
                
                Worksheets(lc_batch).Activate
            End If
            n = n + 1
            Range("L" & n).Activate
        Loop
        
        ' Se da formato a todas las columnas del envelope
        Worksheets(lc_batch).Activate
        Range("K:K").Copy
        Worksheets(envelope_nombre).Activate
        Columns("B:G").Select
        Selection.PasteSpecial Paste:=xlPasteFormats
        'Columns("E:E").Select
        'Selection.PasteSpecial Paste:=xlPasteFormats
        Application.CutCopyMode = False

        ' A las columnas con el número de caso de carga se da formato de número sin decimales
        Columns("C:C").Select
        Selection.NumberFormat = "0"
        Columns("F:F").Select
        Selection.NumberFormat = "0"
        
        ' A las columnas con el nombre del caso de carga se da formato de texto
        Columns("D:D").Select
        Selection.NumberFormat = "@"
        Columns("G:G").Select
        Selection.NumberFormat = "@"
        
        ' Se ajusta el ancho de las columnas
        Columns("A:G").EntireColumn.AutoFit
        Cells(1, 1).Activate
    End If
' -----------------
' FIN DE LA CONSTRUCCIÓN DEL ENVELOPE

    ' Auto guardado
    If Auto_Save = "SI" Then
        ActiveWorkbook.Save
    End If
End Sub




