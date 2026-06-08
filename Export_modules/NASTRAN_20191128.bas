Attribute VB_Name = "NASTRAN_20191128"
Function Integer_to_Nastran(Number As String)
    'Dim tipo As String
    Dim diferencia As Double
    Dim salida As String
    'Dim n As Double
    'Dim m As Double
    Dim longitud As Double
    
    
    ' Caso 0:Celda vacía --> Devuelve 8 espacios en blanco
    If Number = "" Or Number = "-" Or Number = "Blank" Or Number = "BLANK" Or Number = "Empty" Or Number = "EMPTY" Or Number = "Vacio" Or Number = "Vacío" Or Number = "VACIO" Or Number = "VACÍO" Then
        Integer_to_Nastran = "        "
        Exit Function
    End If
    
    
    longitud = Len(Number)
    
    
    ' Caso 1: Números con número de caracteres <=8 --> Se copian tal cual añadiendo espacios al final
    If longitud <= 8 Then
        diferencia = 8 - longitud
        salida = Number
        If diferencia > 0 Then
            For n = 1 To diferencia
                salida = salida & " "
            Next n
        End If
        Integer_to_Nastran = salida
        Exit Function
    End If
    
End Function

Function Real_to_Nastran(Number As String)
    Dim tipo As String
    Dim diferencia As Double
    Dim salida As String
    Dim n As Double
    Dim m As Double
    Dim longitud As Double
    
    ' Caso 0:Celda vacía --> Devuelve 8 espacios en blanco
    If Number = "" Or Number = "-" Or Number = "Blank" Or Number = "BLANK" Or Number = "Empty" Or Number = "EMPTY" Or Number = "Vacio" Or Number = "Vacío" Or Number = "VACIO" Or Number = "VACÍO" Then
        Real_to_Nastran = "        "
        Exit Function
    End If
    
    
    tipo = Entero_o_Real(Number)
    If tipo = "Entero" Then
        Number = Number & "."
    End If
    longitud = Len(Number)

    
    ' Caso 1: Números con número de caracteres <=8 --> Se copian tal cual añadiendo espacios al final
    If longitud <= 8 Then
        diferencia = 8 - longitud
        salida = Number
        If diferencia > 0 Then
            For n = 1 To diferencia
                salida = salida & " "
            Next n
        End If
        Real_to_Nastran = salida
        Exit Function
    End If
    
    
    ' Caso 2: Números con más de 8 caracteres pero son cifras no significativas --> se gana precisión si se hace esto en lugar de aplicar la notación científica
    If Number > 0.1 And Number < 9999999 And tipo = "Real" Then
        Real_to_Nastran = Left(Number, 8)
        Exit Function
    End If
    
    ' Caso 3: Números negativos con más de 8 caracteres pero son cifras no significativas --> se gana precisión si se hace esto en lugar de aplicar la notación científica
    If Number < -0.1 And Number > -999999 And tipo = "Real" Then
        Real_to_Nastran = Left(Number, 8)
        Exit Function
    End If
    
    
    ' Caso 4: Notación científica positiva --> por defecto sale con 8 caracteres
    If longitud > 8 And Number > 0 Then
        Real_to_Nastran = Format(Number, "Scientific")
    End If
    
    ' Caso 5: Notación científica negativa --> por defecto sale con 9 caracteres de modo que le quito la "E" para que se quede en 8 y Nastran se lo traga igualmente
    If longitud > 8 And Number < 0 Then
        Real_to_Nastran = Replace(Format(Number, "Scientific"), "E", "")
    End If

End Function
Function Entero_o_Real(Number As String)
    Dim aux As Double
    
    aux = InStr(Number, ".")
    
    If aux = 0 Then
        Entero_o_Real = "Entero"
    Else
        Entero_o_Real = "Real"
    End If
    
End Function
Sub pruebas_formato()
    Dim Num As String
    Dim aux As String
    
    Num = 123456789
    
    aux = Format(Num, "Scientific")
    
End Sub

Sub Gap_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    'Dim Text_izq As String
    'Dim Text_dch As String
    Dim Node1 As String
    Dim Node2 As String
    Dim i As Double
    Dim N_Pcomp As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    
    'Se graban las variables
    path = ActiveSheet.Range("C2").Value
    file_old = ActiveSheet.Range("C3").Value
    file_mod = ActiveSheet.Range("C4").Value
    N_Pcomp = ActiveSheet.Range("F2").Value
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            Node1 = Mid(textline, 26, 7)
            Node2 = Mid(textline, 34, 7)

            
            For i = 1 To N_Pcomp
                If Node1 = Cells(i + 3, 6).Value Then
                    test = 1
                End If
                If Node2 = Cells(i + 3, 6).Value Then
                    test = 1
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
End Sub

Sub GRID_Comment()
    
    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim GRID_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Nodes As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("C2").Value
    file_old = ActiveSheet.Range("C3").Value
    file_mod = ActiveSheet.Range("C4").Value
    N_Nodes = ActiveSheet.Range("B7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            GRID_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 4)

            
            For i = 1 To N_Nodes
                If GRID_ID = Cells(i + 8, 2).Value And COMIENZO = "GRID" Then
                    test = 1
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Nodes commented"
End Sub
Sub GRID_CD_MOD()
    'AVISO: Esta macro solo funciona en el caso en que los campos 8 y 9 de las tarjetas "GRID" no estén escritos. Cosa que suele ser habitual. Es decir, lo normal es que el último campo escrito sea el "CD" o incluso que este esté en blanco para adoptar el "CD" por defecto que es el "Basic"
    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim GRID_ID As String
    Dim TODO_MENOS_CD As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Nodes As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("I2").Value
    file_old = ActiveSheet.Range("I3").Value
    file_mod = ActiveSheet.Range("I4").Value
    N_Nodes = ActiveSheet.Range("H7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 4)
            If COMIENZO = "GRID" Then 'Solo se recorre la lista de nodos de la hoja Excel si la línea en cuestión es un nodo. Es para aligerar el codigo
            
                GRID_ID = Trim(Mid(textline, 9, 8))
                TODO_MENOS_CD = Left(textline, 48)
                    
                For i = 1 To N_Nodes
                    If GRID_ID = Cells(i + 8, 8).Value Then
                        test = 1
                        Do While Len(TODO_MENOS_CD) < 48 'Este bucle es para rellenar con espacios hasta la posición 48 xq en la 49 empieza el CD
                            TODO_MENOS_CD = TODO_MENOS_CD & " "
                        Loop
                        A.WriteLine (TODO_MENOS_CD & Integer_to_Nastran(Cells(i + 8, 8 + 1).Value))
                        n = n + 1
                    End If
                Next i
            End If
            If test = 0 Then
                A.WriteLine (textline)
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Nodes CD Modified"
End Sub

Sub GRID_CP_MOD()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim GRID_ID As String
    Dim Parte_izq As String
    Dim Parte_dch As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Nodes As Double
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("C2").Value
    file_old = ActiveSheet.Range("C3").Value
    file_mod = ActiveSheet.Range("C4").Value
    N_Nodes = ActiveSheet.Range("B7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 4)
            If COMIENZO = "GRID" Then 'Solo se recorre la lista de nodos de la hoja Excel si la línea en cuestión es un nodo. Es para aligerar el codigo
            
                GRID_ID = Trim(Mid(textline, 9, 8))
                Parte_izq = Left(textline, 16)                      'La cadena de texto que está antes del CP
                Parte_dch = Right(textline, Len(textline) - 24)     'La cadena de texto que está después del CP
                    
                For i = 1 To N_Nodes
                    If GRID_ID = Cells(i + 8, 2).Value Then
                        test = 1
                        A.WriteLine (Parte_izq & Integer_to_Nastran(Cells(i + 8, 2 + 1).Value) & Parte_dch)
                        n = n + 1
                    End If
                Next i
            End If
            If test = 0 Then
                A.WriteLine (textline)
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Nodes CP Modified"
End Sub

Sub CELAS_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim CELAS_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Celas As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("I2").Value
    file_old = ActiveSheet.Range("I3").Value
    file_mod = ActiveSheet.Range("I4").Value
    N_Celas = ActiveSheet.Range("H7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            CELAS_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_Celas
                If CELAS_ID = Cells(i + 8, 8).Value Then
                    If COMIENZO = "CELAS1" Or COMIENZO = "CELAS2" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Celas commented"
End Sub

Sub CBUSH_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim CBUSH_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Cbush As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("O2").Value
    file_old = ActiveSheet.Range("O3").Value
    file_mod = ActiveSheet.Range("O4").Value
    N_Cbush = ActiveSheet.Range("N7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            CBUSH_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 5)

            
            For i = 1 To N_Cbush
                If CBUSH_ID = Cells(i + 8, 14).Value And COMIENZO = "CBUSH" Then
                    test = 1
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Cbush commented"
End Sub

Sub CQUAD4_CTRIA3_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim ELM_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Elm As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("U2").Value
    file_old = ActiveSheet.Range("U3").Value
    file_mod = ActiveSheet.Range("U4").Value
    N_Elm = ActiveSheet.Range("T7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            ELM_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_Elm
                If ELM_ID = Cells(i + 8, 20).Value Then
                    If COMIENZO = "CQUAD4" Or COMIENZO = "CTRIA3" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " CQUAD4 or CTRIA3 elements commented"
End Sub
Sub FORCE_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim FORCE_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Force As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AA2").Value
    file_old = ActiveSheet.Range("AA3").Value
    file_mod = ActiveSheet.Range("AA4").Value
    N_Force = ActiveSheet.Range("Z7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            FORCE_ID = Trim(Mid(textline, 17, 8))
            COMIENZO = Left(textline, 5)

            
            For i = 1 To N_Force
                If FORCE_ID = Cells(i + 8, 26).Value Then
                    If COMIENZO = "FORCE" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Forces commented by Nodes"
End Sub
Sub FORCE_SID_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim FORCE_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Force As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AF2").Value
    file_old = ActiveSheet.Range("AF3").Value
    file_mod = ActiveSheet.Range("AF4").Value
    N_Force = ActiveSheet.Range("AE7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            FORCE_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 5)

            
            For i = 1 To N_Force
                If FORCE_ID = Cells(i + 8, 31).Value Then
                    If COMIENZO = "FORCE" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Forces commented by SID"
End Sub
Sub MOMENT_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim MOMENT_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Moment As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AL2").Value
    file_old = ActiveSheet.Range("AL3").Value
    file_mod = ActiveSheet.Range("AL4").Value
    N_Moment = ActiveSheet.Range("AK7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            MOMENT_ID = Trim(Mid(textline, 17, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_Moment
                If MOMENT_ID = Cells(i + 8, 37).Value Then
                    If COMIENZO = "MOMENT" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Moments commented by Nodes"
End Sub
Sub MOMENT_SID_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim MOMENT_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Moment As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AQ2").Value
    file_old = ActiveSheet.Range("AQ3").Value
    file_mod = ActiveSheet.Range("AQ4").Value
    N_Moment = ActiveSheet.Range("AP7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            MOMENT_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_Moment
                If MOMENT_ID = Cells(i + 8, 42).Value Then
                    If COMIENZO = "MOMENT" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Moments commented by SID"
End Sub

Sub PLOAD4_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim PLOAD4_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_PLOAD4 As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AW2").Value
    file_old = ActiveSheet.Range("AW3").Value
    file_mod = ActiveSheet.Range("AW4").Value
    N_PLOAD4 = ActiveSheet.Range("AV7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            PLOAD4_ID = Trim(Mid(textline, 17, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_PLOAD4
                If PLOAD4_ID = Cells(i + 8, 48).Value Then
                    If COMIENZO = "PLOAD4" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " PLOAD4 commented by Elements"
End Sub

Sub PLOAD4_SID_Comment()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim PLOAD4_ID As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_PLOAD4 As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("BB2").Value
    file_old = ActiveSheet.Range("BB3").Value
    file_mod = ActiveSheet.Range("BB4").Value
    N_PLOAD4 = ActiveSheet.Range("BA7").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            PLOAD4_ID = Trim(Mid(textline, 9, 8))
            COMIENZO = Left(textline, 6)

            
            For i = 1 To N_PLOAD4
                If PLOAD4_ID = Cells(i + 8, 53).Value Then
                    If COMIENZO = "PLOAD4" Then
                        test = 1
                    End If
                End If
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            Else
                A.WriteLine ("$" & textline)
                n = n + 1
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " PLOAD4 commented by SID"
End Sub

Sub pcomp_modifier()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim Text_izq As String
    Dim Text_dch As String
    Dim i As Double
    Dim N_Pcomp As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    
    'Se graban las variables
    path = ActiveSheet.Range("C2").Value
    file_old = ActiveSheet.Range("C3").Value
    file_mod = ActiveSheet.Range("C4").Value
    N_Pcomp = ActiveSheet.Range("F2").Value
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
        
            Line Input #1, textline
            Text_izq = Left(textline, 3 * 8)
            If Len(textline) > 3 * 8 Then
                Text_dch = Right(textline, Len(textline) - 3 * 8)
            End If
            
            
            For i = 1 To N_Pcomp
                Left_string = "CQUAD4   " & Cells(i + 3, 6).Value & " " & Cells(i + 3, 7).Value
                New_Left_string = "CQUAD4   " & Cells(i + 3, 6).Value & " " & Cells(i + 3, 8).Value
                
                If Text_izq = Left_string Then
                    A.WriteLine (New_Left_string & Text_dch)
                    test = 1
                End If
                
            Next i
            
            If test = 0 Then
                A.WriteLine (textline)
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
End Sub

Sub grids_to_file()
    'Se declaran las variables
    Dim path As String
    Dim file As String
    Dim N_Grids As Double
    path = ActiveSheet.Range("I2").Value
    file = ActiveSheet.Range("I3").Value
    N_Grids = ActiveSheet.Range("I6").Value
    'Se crea el fichero modificado
    NewFile = (path & "\" & file)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)

    Dim i As Double
    For i = 1 To N_Grids
        A.WriteLine (Cells(i + 6, 9).Value)
    Next i
    'Se cierra el fichero
    A.Close
    
End Sub

Sub Nastran_bdf_Split()
    
    'Split any Column into 8 character cells(previous duplication of the column)-->Keeps un-splitted values
    Range("A:A").Select
    colu = ActiveCell.column 'en que columna estas
    cell_row = Selection.Row 'en que fila estas
    Selection.Copy 'copias columna
    
    Columns(colu + 1).Select 'seleccionas siguiente columna
    ActiveSheet.Paste 'pegas

    'separar en columnas
    Selection.TextToColumns Destination:=Cells(cell_row, colu + 1), DataType:=xlFixedWidth, _
        FieldInfo:=Array(Array(0, 1), Array(8, 1), Array(16, 1), Array(24, 1), Array(32, 1), _
        Array(40, 1), Array(48, 1), Array(56, 1), Array(64, 1), Array(72, 1), Array(80, 1)), _
        TrailingMinusNumbers:=True
    
    Columns("B:K").Select
    Columns("B:K").EntireColumn.AutoFit
    Range("A1").Select

End Sub

Function Nastran_Card_Maker(fields As Range)

Dim i As Long
Dim n As Long
Dim v_fields() As String

n = maxrange(fields)
ReDim v_fields(1 To n)

For i = 1 To n
    v_fields(i) = fields(i)
Next i

Nastran_Card_Maker = ""
For i = 1 To n
    Nastran_Card_Maker = Nastran_Card_Maker & ExcelSpcsNastran(v_fields(i))
Next i

End Function

Sub Multiple_Search()

'Se declaran las variables
Dim path As String
Dim file As String
Dim v_find() As String
Dim v_discard() As String
Dim n_find As Long
Dim n_discard As Long
Dim linea_num As Double ' para saber en que linea encuentra lo que buscas
Dim finds_counter As Long ' cuenta cuantos de los criterios encuentra
Dim discards_counter As Long ' igual con los que quieres que no encuentre
Dim cont_aciertos As Long

'eliminar resultados previos
Range("J2:K2").Select
Range(Selection, Selection.End(xlDown)).Select
'Range(Selection, Selection.End(xlDown)).Select 'si en vez de llegar hasta rango relleno se quiere llegar a final de la hoja, duplicar este comando2
Selection.ClearContents

'Se graban las variables
path = ActiveSheet.Range("C2").Value
file = ActiveSheet.Range("C3").Value
n_find = ActiveSheet.Range("F1").Value
n_discard = ActiveSheet.Range("H1").Value



'Criterios de busqueda
ReDim v_find(1 To n_find)
ReDim v_discard(1 To n_discard)

For i = 1 To n_find
    v_find(i) = Cells(1 + i, 6).Value
Next i

For i = 1 To n_discard
    v_discard(i) = Cells(1 + i, 8).Value
Next i



linea_num = 0 'Inicializar contador de linea
cont_aciertos = 0 'cuantas lineas coinciden con todos los criterios


'Se abre el fichero en el que se busca
Open (path & "\" & file) For Input As #1

While Not EOF(1)

    Line Input #1, textline
    linea_num = linea_num + 1
    
    finds_counter = 0
    For i = 1 To UBound(v_find)
        
        If InStr(1, textline, v_find(i), vbTextCompare) > 0 Then '>0 porque da la posicion en la que lo ecuentra
            
            finds_counter = finds_counter + 1
        
        End If
    
    Next i
        
    'En caso de que encuentre todos los criterios, te vale si no estan las cosas que no quieres que esten
    discards_counter = 0
    If finds_counter = UBound(v_find) Then
        
        For i = 1 To UBound(v_discard)
            
            If InStr(1, textline, v_discard(i), vbTextCompare) > 0 Then '=0 significa no lo encuentra

                discards_counter = discards_counter + 1
            End If
        Next i
        
    End If
    
    
    
    If finds_counter = UBound(v_find) And discards_counter = 0 Then
        
        ActiveSheet.Cells(2 + cont_aciertos, 10) = textline
        ActiveSheet.Cells(2 + cont_aciertos, 11) = linea_num
        cont_aciertos = cont_aciertos + 1
    
    End If

Wend
    
'Se cierra el fichero fuente
Close #1


End Sub

Function split_by_length_horiz(text As String, length As Long)

Dim long_text As Long
Dim resto As Long
Dim v_trozos() As String



long_text = Len(text)
If long_text <= length Then

    ReDim v_trozos(1 To 1)
    v_trozos(1) = text
    
Else
    resto = long_text Mod length
    
    num_trozos = WorksheetFunction.RoundDown(long_text / length, 0)
    ReDim v_trozos(1 To num_trozos + 1)
    
    For i = 1 To (num_trozos + 1) * length Step length 'num de trozos redondeados +1 del resto que falta
        
        If i = 1 Then ' primer cacho
        
            v_trozos(i) = Mid(text, i, length)

        ElseIf i = (num_trozos) * length + 1 Then  'ultimo cacho
        
            v_trozos(((i - 1) / length) + 1) = Mid(text, i, resto)
            
        Else 'resto de cachos
        
            v_trozos(((i - 1) / length) + 1) = Mid(text, i, length)
            
        End If
        
    Next i

End If

split_by_length_horiz = v_trozos

End Function

Function split_by_length_vert(text As String, length As Long)

Dim long_text As Long
Dim resto As Long
Dim v_trozos() As String



long_text = Len(text)
If long_text <= length Then

    ReDim v_trozos(1 To 1)
    v_trozos(1) = text
    
Else
    resto = long_text Mod length
    
    num_trozos = WorksheetFunction.RoundDown(long_text / length, 0)
    ReDim v_trozos(1 To num_trozos + 1)
    
    For i = 1 To (num_trozos + 1) * length Step length 'num de trozos redondeados +1 del resto que falta
        
        If i = 1 Then ' primer cacho
        
            v_trozos(i) = Mid(text, i, length)

        ElseIf i = (num_trozos) * length + 1 Then  'ultimo cacho
        
            v_trozos(((i - 1) / length) + 1) = Mid(text, i, resto)
            
        Else 'resto de cachos
        
            v_trozos(((i - 1) / length) + 1) = Mid(text, i, length)
            
        End If
        
    Next i

End If

split_by_length_vert = WorksheetFunction.Transpose(v_trozos)

End Function

Function Concat_text(parts As Range, delimiter As String)

Dim i As Long
Dim n As Long
Dim v_parts() As String

n = maxrange(parts)
ReDim v_parts(1 To n)

'almacenar rango en vector
For i = 1 To n
    v_parts(i) = parts(i)
Next i


Concat_text = ""
For i = 1 To n - 1
    Concat_text = Concat_text & v_parts(i) & delimiter
Next i

    Concat_text = Concat_text & v_parts(n)


End Function

Sub ELEMENT_PID_MOD()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim ELEMENT_ID As String
    Dim Parte_izq As String
    Dim Parte_dch As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Elements As Double
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("Q2").Value
    file_old = ActiveSheet.Range("Q3").Value
    file_mod = ActiveSheet.Range("Q4").Value
    N_Elements = ActiveSheet.Range("P7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 6)
            If COMIENZO = "CQUAD4" Or COMIENZO = "CTRIA3" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el codigo
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))              'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                Parte_izq = Left(textline, 16)                      'La cadena de texto que está antes del PID. Hay dos campos antes del CP --> 2·8=16 caracteres
                Parte_dch = Right(textline, Len(textline) - 24)     'La cadena de texto que está después del PID. Se quedan 3 campos después del CP --> 3·8=24 caracteres, es decir nos quedamos con todo menos los 24 prieros caracteres
                    
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 16).Value Then     'La lista de elementos está en la columna 16 de la hoja Excel
                        test = 1
                        A.WriteLine (Parte_izq & Integer_to_Nastran(Cells(i + 8, 16 + 1).Value) & Parte_dch)
                        n = n + 1
                    End If
                Next i
            End If
            If test = 0 Then
                A.WriteLine (textline)
            End If
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Elements PID Modified"
End Sub

Sub ELEMENT_OFFSET_MOD()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim ELEMENT_ID As String
    Dim TODO_MENOS_OFFSET As String
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Elements As Double
    Dim Left_string As String
    Dim New_Left_string As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("W2").Value
    file_old = ActiveSheet.Range("W3").Value
    file_mod = ActiveSheet.Range("W4").Value
    N_Elements = ActiveSheet.Range("V7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 6)
            
            '--------
            If COMIENZO = "CQUAD4" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el codigo
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))              'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                TODO_MENOS_OFFSET = Left(textline, 64)              'Nos quedamos con los 8 campos anteriores al OFFSET --> 8·8=64 caracteres
                    
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 22).Value Then     'La lista de elementos está en la columna 22 de la hoja Excel
                        test = 1
                        Do While Len(TODO_MENOS_OFFSET) < 64 'Este bucle es para rellenar con espacios hasta la posición 64 xq en la 65 empieza el OFFSET
                            TODO_MENOS_OFFSET = TODO_MENOS_OFFSET & " "
                        Loop
                        A.WriteLine (TODO_MENOS_OFFSET & Real_to_Nastran(Cells(i + 8, 22 + 1).Value))   'Ojo que el OFFSET es un campo real asique esta función llama a la función  "Real_to_Nastran", no a la "Integer_to_Nastran"
                        n = n + 1
                    End If
                Next i
            End If
            '--------
            
            '--------
            If COMIENZO = "CTRIA3" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento CTRIA3. Es para aligerar el codigo. El motivo de diferenciarlo entre CQUAD4 y CTRIA3 es que al tener un nodo menos los campos están desplazados hacia la izquierda
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))              'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                TODO_MENOS_OFFSET = Left(textline, 56)              'Nos quedamos con los 7 campos anteriores al OFFSET --> 7·8=56 caracteres
                    
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 22).Value Then     'La lista de elementos está en la columna 22 de la hoja Excel
                        test = 1
                        Do While Len(TODO_MENOS_OFFSET) < 56 'Este bucle es para rellenar con espacios hasta la posición 56 xq en la 57 empieza el OFFSET
                            TODO_MENOS_OFFSET = TODO_MENOS_OFFSET & " "
                        Loop
                        A.WriteLine (TODO_MENOS_OFFSET & Real_to_Nastran(Cells(i + 8, 22 + 1).Value))   'Ojo que el OFFSET es un campo real asique esta función llama a la función  "Real_to_Nastran", no a la "Integer_to_Nastran"
                        n = n + 1
                    End If
                Next i
            End If
            '--------
            
            '--------
            If test = 0 Then            'Si la línea no es ni un CQUAD4 ni un CTRIA3 entonces se escribe sin ningún cambio
                A.WriteLine (textline)
            End If
            '--------
            
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Elements OFFSET Modified"
End Sub

Sub ELEMENT_MCID_THETA_MOD()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim ELEMENT_ID As String
    Dim Parte_izq As String
    Dim Parte_dch As String
    Dim Total_Len As Integer
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Elements As Double
    Dim Field_type As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AC2").Value
    file_old = ActiveSheet.Range("AC3").Value
    file_mod = ActiveSheet.Range("AC4").Value
    N_Elements = ActiveSheet.Range("AB7").Value
    Field_type = ActiveSheet.Range("AC8").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 6)
            
            '--------
            If COMIENZO = "CQUAD4" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el código
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))                  'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                Parte_izq = Left(textline, 56)                          'La cadena de texto que está antes del MCID/THETA. Hay 7 campos antes del MCID/THETA --> 7·8=56 caracteres
                Total_Len = Len(textline)                               'La longitud total de la cadena de caracteres para saber si hay o no hay "parte derecha"
                If Total_Len > 64 Then
                    Parte_dch = Right(textline, Len(textline) - 64)     'La cadena de texto que está después del MCID/THETA. Se quedan 1 campo después del MCID/THETA --> 1·8=8 caracteres, es decir nos quedamos con todo menos los 8·8=64 prieros caracteres incluyendo el propio campo de MCID/THETA. Pero solamente si de verdad existen. De ahí la existencia del condicional
                End If
                
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 28).Value Then         'La lista de elementos está en la columna 28 de la hoja Excel
                        test = 1
                        Do While Len(Parte_izq) < 56 'Este bucle es para rellenar con espacios hasta la posición 56 xq en la 57 empieza el OFFSET
                            Parte_izq = Parte_izq & " "
                        Loop
                        If Field_type = "NUEVO MCID" Then
                            A.WriteLine (Parte_izq & Integer_to_Nastran(Cells(i + 8, 28 + 1).Value) & Parte_dch)
                        End If
                        If Field_type = "NUEVO THETA" Then
                            A.WriteLine (Parte_izq & Real_to_Nastran(Cells(i + 8, 28 + 1).Value) & Parte_dch)
                        End If
                        n = n + 1
                    End If
                Next i
            End If
            '--------
            

            '--------
            If COMIENZO = "CTRIA3" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el código
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))                  'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                Parte_izq = Left(textline, 48)                          'La cadena de texto que está antes del MCID/THETA. Hay 6 campos antes del MCID/THETA --> 6·8=48 caracteres
                Total_Len = Len(textline)                               'La longitud total de la cadena de caracteres para saber si hay o no hay "parte derecha"
                If Total_Len > 56 Then
                    Parte_dch = Right(textline, Len(textline) - 56)     'La cadena de texto que está después del MCID/THETA. Se quedan 1 campo después del MCID/THETA --> 1·8=8 caracteres, es decir nos quedamos con todo menos los 7·8=56 primeros caracteres incluyendo el propio campo de MCID/THETA. Pero solamente si de verdad existen. De ahí la existencia del condicional
                End If
                
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 28).Value Then         'La lista de elementos está en la columna 28 de la hoja Excel
                        test = 1
                        Do While Len(Parte_izq) < 48 'Este bucle es para rellenar con espacios hasta la posición 48 xq en la 49 empieza el OFFSET
                            Parte_izq = Parte_izq & " "
                        Loop
                        If Field_type = "NUEVO MCID" Then
                            A.WriteLine (Parte_izq & Integer_to_Nastran(Cells(i + 8, 28 + 1).Value) & Parte_dch)
                        End If
                        If Field_type = "NUEVO THETA" Then
                            A.WriteLine (Parte_izq & Real_to_Nastran(Cells(i + 8, 28 + 1).Value) & Parte_dch)
                        End If
                        n = n + 1
                    End If
                Next i
            End If
            '--------
            
            
            '--------
            If test = 0 Then
                A.WriteLine (textline)
            End If
            '--------
            
            Parte_dch = ""
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Elements MCID/THETA Modified"
End Sub

Sub ELEMENT_INVERT_NORMAL()

    'Se declaran las variables
    Dim path As String
    Dim file_old As String
    Dim file_mod As String
    Dim ELEMENT_ID As String
    Dim Parte_izq As String
    Dim Parte_dch As String
    Dim Total_Len As Integer
    Dim COMIENZO As String
    Dim i As Double
    Dim N_Elements As Double
    Dim Nodo_1 As String
    Dim Nodo_2 As String
    Dim Nodo_3 As String
    Dim Nodo_4 As String
    Dim test As Integer
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("AI2").Value
    file_old = ActiveSheet.Range("AI3").Value
    file_mod = ActiveSheet.Range("AI4").Value
    N_Elements = ActiveSheet.Range("AH7").Value

    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"
    
    
    'Se abre el fichero fuente
    Open (path & "\" & file_old) For Input As #1

    'Se crea el fichero modificado
    NewFile = (path & "\" & file_mod)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
        While Not EOF(1)
        
            Line Input #1, textline
            COMIENZO = Left(textline, 6)
            
            '--------
            If COMIENZO = "CQUAD4" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el código
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))                  'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                Parte_izq = Left(textline, 24)                          'La cadena de texto que está antes del primer nodo. Hay 3 campos antes del primer nodo --> 3·8=24 caracteres
                Total_Len = Len(textline)                               'La longitud total de la cadena de caracteres para saber si hay o no hay "parte derecha"
                If Total_Len > 56 Then
                    Parte_dch = Right(textline, Len(textline) - 56)     'La cadena de texto que está después del ultimo nodo. Podría haber 1 ó 2 campos después del último nodo --> Es decir nos quedamos con todo menos los 7·8=56 prieros caracteres incluyendo el último nodo. Pero solamente si de verdad existen. De ahí la existencia del condicional
                End If
                
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 34).Value Then         'La lista de elementos está en la columna 34 de la hoja Excel
                        test = 1
                        Nodo_1 = Trim(Mid(textline, 25, 8))
                        Nodo_2 = Trim(Mid(textline, 33, 8))
                        Nodo_3 = Trim(Mid(textline, 41, 8))
                        Nodo_4 = Trim(Mid(textline, 49, 8))
                        A.WriteLine (Parte_izq & Integer_to_Nastran(Nodo_4) & Integer_to_Nastran(Nodo_3) & Integer_to_Nastran(Nodo_2) & Integer_to_Nastran(Nodo_1) & Parte_dch)
                        n = n + 1
                    End If
                Next i
            End If
            '--------


            '--------
            If COMIENZO = "CTRIA3" Then 'Solo se recorre la lista de elementos de la hoja Excel si la línea en cuestión es un elemento. Es para aligerar el código
            
                ELEMENT_ID = Trim(Mid(textline, 9, 8))                  'Nos quedamos con los 8 caractereres seguidos al noveno que corresponden al ID del elemento
                Parte_izq = Left(textline, 24)                          'La cadena de texto que está antes del primer nodo. Hay 3 campos antes del primer nodo --> 3·8=24 caracteres
                Total_Len = Len(textline)                               'La longitud total de la cadena de caracteres para saber si hay o no hay "parte derecha"
                If Total_Len > 48 Then
                    Parte_dch = Right(textline, Len(textline) - 48)     'La cadena de texto que está después del ultimo nodo. Podría haber 1 ó 2 campos después del último nodo --> Es decir nos quedamos con todo menos los 6·8=48 prieros caracteres incluyendo el último nodo. Pero solamente si de verdad existen. De ahí la existencia del condicional
                End If
                
                For i = 1 To N_Elements
                    If ELEMENT_ID = Cells(i + 8, 34).Value Then         'La lista de elementos está en la columna 34 de la hoja Excel
                        test = 1
                        Nodo_1 = Trim(Mid(textline, 25, 8))
                        Nodo_2 = Trim(Mid(textline, 33, 8))
                        Nodo_3 = Trim(Mid(textline, 41, 8))
                        'Nodo_4 = Trim(Mid(textline, 49, 8))
                        A.WriteLine (Parte_izq & Integer_to_Nastran(Nodo_3) & Integer_to_Nastran(Nodo_2) & Integer_to_Nastran(Nodo_1) & Parte_dch)
                        n = n + 1
                    End If
                Next i
            End If
            '--------

            
            '--------
            If test = 0 Then
                A.WriteLine (textline)
            End If
            '--------
            
            Parte_dch = ""
            test = 0
        Wend
        
    'Se cierra el fichero fuente
    Close #1
    
    'Se cierra el fichero modificado
    A.Close
    
    Application.StatusBar = n & " Elements inverted normals"
End Sub

Sub Linear_Gap_Creation()

    'Se declaran las variables
    Dim path As String
    Dim file_name As String
    Dim N_Contacts As Double
    Dim i As Double
    Dim j As Double
    Dim n As Double
    
    'Se graban las variables
    path = ActiveSheet.Range("U2").Value
    file_name = ActiveSheet.Range("U3").Value
    N_Contacts = ActiveSheet.Range("U4").Value
    
    Application.DisplayStatusBar = True
    Application.StatusBar = "Running"

    'Se crea el fichero
    NewFile = (path & "\" & file_name)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)
    
    n = 0
    
    For i = 1 To N_Contacts
        For j = 1 To 6
            A.WriteLine (Cells(i + 3, j + 12).Value)
        Next j
        n = n + 1
    Next i
    
    'Se cierra el fichero
    A.Close
    
    Application.StatusBar = n & " Linear Gaps created"
End Sub


Function coincidir_nodos(id1 As Range, X1 As Range, Y1 As Range, Z1 As Range, id2 As Range, X2 As Range, Y2 As Range, Z2 As Range)

    'Find the length of the range 1
    Dim l1 As Integer
    l1 = maxrange(id1)

    'Find the length of the range 2
    Dim l2 As Integer
    l2 = maxrange(id2)

    'New vector with le length of the range
    'First column (0) id of the closest node. Second column (1) distance
    ReDim vec(l2, 1) As Double

    
    Dim i As Double
    Dim j As Double
    Dim distance As Double
    Dim distance_temp As Double
    Dim new_node As Double
    
    For i = 1 To l2
        distance = 1E+300                           'Se reestablece la variable empezando en un valor absurdamente elevado
        distance_temp = 1E+300                      'Esto no sería necesario, pero simplifica la trazabilidad
        For j = 1 To l1
            distance_temp = ((X2(i) - X1(j)) ^ 2 + (Y2(i) - Y1(j)) ^ 2 + (Z2(i) - Z1(j)) ^ 2) ^ (0.5)   'Se calcula la distancia
            If distance_temp < distance Then        'Si la nueva distancia es menor que la mínima almacenada, la reemplazados y almacenamos el valor de ese nodo
                distance = distance_temp
                new_node = id1(j)
            End If
        Next j
        vec(i - 1, 0) = new_node                    'Sacamos el valor del nuevo nodo
        vec(i - 1, 1) = distance                    'Sacamos el valor de la distancia que en principio, serviría para chequeo
    Next i
        
    coincidir_nodos = vec                           'Se recupera el vector de resultado
End Function
