Attribute VB_Name = "Patran_20190723"
Function Compactador_Patran(ran As Range)
    Dim Num_ini As Double
    Dim Num_fin As Double
    Dim contador As Double
    Dim salida As String
    Dim i As Double
    Dim Num As Double
    Num = UBound(ran(), 1)

    For i = 1 To Num
    
        If i = 1 Then
            Num_ini = ran(1, 1)
        End If
        
        If ran(i + 1, 1) = ran(i, 1) + 1 Then
            Num_fin = ran(i + 1, 1)
            contador = contador + 1
        End If
        
        If ran(i + 1, 1) <> ran(i, 1) + 1 And contador = 0 Then
            salida = salida & " " & Num_ini
            Num_ini = ran(i + 1, 1)
        End If
        
        If ran(i + 1, 1) <> ran(i, 1) + 1 And contador <> 0 Then
            salida = salida & " " & Num_ini & ":" & Num_fin
            Num_ini = ran(i + 1, 1)
            contador = 0
        End If
        
    Next i
    
    Compactador_Patran = salida

End Function

Sub Patran_List()
    'OBSOLETE --> SEE "Patran_List_Complex" SUBROUTINE
    'This subroutine allows to create a list of elements starting from a Patran simple list
    
    'Reading of the list
    Dim list As String
    list = Range("C2").Value
    
    'Find the position of the separator of the list
    Dim pos_sep As Double
    pos_sep = InStr(1, list, ":")
    
    'Separation of the list
    Dim start As Double
    Dim ending As Double
    start = Left(list, pos_sep - 1)
    ending = Right(list, Len(list) - pos_sep)
    
    'Write the new list
    Range("C4").Activate
    Dim i As Double
    Dim n As Double
    For i = start To ending
        'Cells(i + 3, 3).Activate
        Cells(n + 4, 3) = i
        n = n + 1
    Next i
    
End Sub

Sub CLEAR_PATRAN_LIST()
    'Clears the result list of the PATRAN_LIST
    Range("C4").Select
    Range(Selection, Selection.End(xlDown)).Select
    Selection.ClearContents
    Range("A1").Select
End Sub

Sub Patran_List_Complex()
    'This subroutine allows to create a list of elements starting from a Patran complex list
    
    Application.ScreenUpdating = False              'Disable screen update for speed up
    
    Dim flag As Integer
    flag = 0
    If Application.Calculation = xlAutomatic Then   'If the book is in automatic calculation mode then
        flag = 1                                    'The current mode is registered for later restore
        Application.Calculation = xlManual          'The calculation is turned into manual mode for speed up
    End If
    
    Dim list As String
    list = Range("C2").Value                        'Reading of the list
    
    list_split = Split(list)                        'Split of the list
    
    Dim N_lists As Double                           'Find the number of sub-lists in the list
    N_lists = UBound(list_split, 1) '+ 1
    
    Range("C4").Activate
    
    Dim m As Double
    Dim i As Double
    Dim n As Double
    Dim k As Double
    Dim kk As Double
    Dim start As Double
    Dim ending As Double
    Dim stepp As Double
    Dim N_Collumn As Double 'Find if there is one or two ":"
    Dim sub_list As String
    Dim sub_sub_list As String
    Dim pos_sep As Double
    Dim long_paso As Double
    
    
    'LOOP
    '--------------------------------------------------------------
    For m = 0 To N_lists
        sub_list = list_split(m)
        N_Collumn = Len(sub_list) - Len(Replace(sub_list, ":", ""))
        If N_Collumn = 0 Then 'Se trata de un elemento único, no de una lista
            Cells(n + 4, 3) = sub_list
            n = n + 1
            k = k + 1                                                   'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
            If k = 1000 Then                                            'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                kk = kk + 1000                                          'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                Application.StatusBar = kk & " Elements listed"         'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                k = 0                                                   'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
            End If                                                      'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
        End If
        If N_Collumn > 0 Then 'Se trata de una lista, con o sin paso intermedio
            If N_Collumn = 1 Then 'Se trata de una lista sin paso intermedio (paso=1)
                temp = Split(sub_list, ":")
                 start = temp(0)
                 ending = temp(1)
                 stepp = 1
            ElseIf N_Collumn = 2 Then 'Se trata de una lista con paso intermedio distinto de "1"
                temp = Split(sub_list, ":") 'calculamos la longitud de caracteres del paso: 0-9 --> 1 caracter; 10-99 --> 2 caracteres; 100-999 --> 3 caracteres; ...
                long_paso = Len(temp(2))    'calculamos la longitud de caracteres del paso: 0-9 --> 1 caracter; 10-99 --> 2 caracteres; 100-999 --> 3 caracteres; ...
                start = temp(0)
                ending = temp(1)
                stepp = temp(2)
            End If

            If stepp > 0 Then
                For i = start To ending Step stepp
                    Cells(n + 4, 3) = i
                    n = n + 1
                    k = k + 1                                           'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                    If k = 1000 Then                                    'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        kk = kk + 1000                                  'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        Application.StatusBar = kk & " Elements listed" 'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        k = 0                                           'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                    End If                                              'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                Next i
            End If
            
            If stepp < 0 Then
                For i = start To ending Step stepp
                    Cells(n + 4, 3) = i
                    n = n + 1
                    k = k + 1                                           'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                    If k = 1000 Then                                    'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        kk = kk + 1000                                  'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        Application.StatusBar = kk & " Elements listed" 'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                        k = 0                                           'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                    End If                                              'Este bloque se repite 3 veces. Solo sirve para ir actualizando el medidor de actividad por si la lista es muy larga. Si se borra este bloque las 3 veces la macro funciona igual
                Next i
            End If
            
        End If
    Next m
    '--------------------------------------------------------------
    'END LOOP
    
    
    If flag = 1 Then                                'The book was in automatic calculation mode
        Application.Calculation = xlAutomatic       'So automatic calculation mode is restored
    End If
    
    Application.ScreenUpdating = True               'Enable screen update
    Application.DisplayStatusBar = True             'Enable information bar
    Application.StatusBar = n & " Elements listed"  'Report the number of elements listed in the information bar
    
End Sub
