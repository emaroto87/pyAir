Attribute VB_Name = "BD"
'Function maxrange(vec As Range)
'    ' Returns the max range of an array MAX(rows;cols)
'    ' Declaration of variables
'    Dim L As Double
'    Dim rows As Double
'    Dim cols As Double
'
'    On Error GoTo check
'
'    ' Search for the number or rows and collunms in the array
'    rows = UBound(vec(), 1)
'    cols = UBound(vec(), 2)
'
'
'    ' Get the max between rows and collunms
'    If cols > rows Then
'        L = cols
'    Else
'        L = rows
'    End If
'
'
'' If range is 1x1 UBound returns an error so result is "1" and exit function
'check:
'    If rows = Empty Then
'        maxrange = 1
'        Exit Function
'    End If
'' Else, return the max between rows and cols
'
'
'    maxrange = L
'End Function

Function maxbd(mat As Range, campo As String, campos As Range, valores As Range)
    ' This function imitales DMAX native Excel function
    Dim n As Double
    Dim m As Double
    Dim j As Double
            
    
    ' Find the number of Fields
    Dim ncampos As Double
    ncampos = maxrange(campos)
    
    ' Find the number of rows and columns of the matrix
    Dim rows As Double
    Dim cols As Double
    rows = UBound(mat(), 1)
    cols = UBound(mat(), 2)

    ' Search for the column with the result
    Dim column As Double
    For n = 1 To cols
        If mat(1, n) = campo Then
            column = n
        End If
    Next n

    ' Start MAX with the minimun Excel value
    Dim MAX As Double
    MAX = -1.79769313486231E+308


    ' Start an array containing the position of the fields
    ReDim posicioncampos(ncampos)
    For m = 0 To (ncampos - 1)
        For n = 1 To cols
            If mat(1, n) = campos(m + 1) Then
                posicioncampos(m) = n
            End If
        Next n
    Next m

    ' Find the max
    For n = 2 To rows
        For m = 0 To (ncampos - 1)
            If mat(n, posicioncampos(m)) = valores(m + 1) Then
                j = j + 1
            End If
        Next m
        If j = ncampos Then
            If mat(n, column) <> Empty Then
                If mat(n, column) > MAX Then
                    MAX = mat(n, column)
                End If
            End If
        End If
        j = 0
    Next n
        
        
    ' If no match available, an #N/A error is returned
    If MAX = -1.79769313486231E+308 Then
        maxbd = CVErr(xlErrNA)
        Exit Function
    End If
        
        
    maxbd = MAX
End Function


Function minbd(mat As Range, campo As String, campos As Range, valores As Range)
    ' This function imitales DMin native Excel function
    Dim n As Double
    Dim m As Double
    Dim j As Double
            
    
    ' Find the number of Fields
    Dim ncampos As Double
    ncampos = maxrange(campos)
    
    ' Find the number of rows and columns of the matrix
    Dim rows As Double
    Dim cols As Double
    rows = UBound(mat(), 1)
    cols = UBound(mat(), 2)

    ' Search for the column with the result
    Dim column As Double
    For n = 1 To cols
        If mat(1, n) = campo Then
            column = n
        End If
    Next n

    ' Start MIN with the maximun Excel value
    Dim MIN As Double
    MIN = 1.79769313486231E+308


    ' Start an array containing the position of the fields
    ReDim posicioncampos(ncampos)
    For m = 0 To (ncampos - 1)
        For n = 1 To cols
            If mat(1, n) = campos(m + 1) Then
                posicioncampos(m) = n
            End If
        Next n
    Next m

    ' Find the min
    For n = 2 To rows
        For m = 0 To (ncampos - 1)
            If mat(n, posicioncampos(m)) = valores(m + 1) Then
                j = j + 1
            End If
        Next m
        If j = ncampos Then
            If mat(n, column) <> Empty Then
                If mat(n, column) < MIN Then
                    MIN = mat(n, column)
                End If
            End If
        End If
        j = 0
    Next n
        
        
    ' If no match available, an #N/A error is returned
    If MIN = 1.79769313486231E+308 Then
        minbd = CVErr(xlErrNA)
        Exit Function
    End If
        
        
    minbd = MIN
End Function



