Attribute VB_Name = "Brackets"
Sub Genera_Bracket_T()

    '-----------------------------------------------------
    'Variable declaration
    Dim i As Double
    Dim j As Double
    Dim Altura As Double
    Dim Diametro As Double
    Dim Young As Double
    Dim Poisson As Double
    Dim Densidad As Double
    Dim Masa1 As String
    Dim X1 As Double
    Dim Y1 As Double
    Dim M1 As Double
    Dim Masa2 As String
    Dim X2 As Double
    Dim Y2 As Double
    Dim M2 As Double
    Dim Ruta As String
    Dim Nombre_fichero As String
    '-----------------------------------------------------
    
    
    '-----------------------------------------------------
    'Variable reading
    Ruta = ActiveSheet.Range("C2")
    Nombre_fichero = ActiveSheet.Range("C3")
    
    Altura = ActiveSheet.Range("C5")
    Diametro = ActiveSheet.Range("C6")
    
    Young = ActiveSheet.Range("C8")
    Poisson = ActiveSheet.Range("C9")
    Densidad = ActiveSheet.Range("C10")
    
    Masa1 = ActiveSheet.Range("C12")
    X1 = ActiveSheet.Range("C13")
    Y1 = ActiveSheet.Range("C14")
    M1 = ActiveSheet.Range("C15")
    
    Masa2 = ActiveSheet.Range("C17")
    X2 = ActiveSheet.Range("C18")
    Y2 = ActiveSheet.Range("C19")
    M2 = ActiveSheet.Range("C20")
    '-----------------------------------------------------
    
    
    '-----------------------------------------------------
    'Check for already existing file
    Dim test As Boolean
    test = False
    If existe_archivo(Ruta & "\" & Nombre_fichero & ".bdf") = True Then
        i = 1
        Do While test = False
            file_name_temp = file_name & i
            If existe_archivo(path & "\" & Nombre_fichero & i & ".bdf") = False Then
                file_name = file_name_temp
                test = True
            End If
            i = i + 1
        Loop
    End If
    '-----------------------------------------------------
    

    '-----------------------------------------------------
    'Create .bdf archive
    NewFile = (Ruta & "\" & Nombre_fichero & ".bdf")
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set A = FSO.createtextfile(NewFile, True)

    'File Header
    A.WriteLine ("$####################################################################")
    A.WriteLine ("$#                PARAMETRIC T BRACKET MODEL CREATOR                #")
    A.WriteLine ("$####################################################################")
    A.WriteLine ("")
    A.WriteLine ("SOL 103")
    A.WriteLine ("CEND")
    A.WriteLine ("")
    A.WriteLine ("TITLE = MSC.Nastran job created on " & Now)
    A.WriteLine ("ECHO = NONE")
    A.WriteLine ("SUBCASE 1")
    A.WriteLine ("   SUBTITLE=Default")
    A.WriteLine ("   METHOD = 1")
    A.WriteLine ("   SPC = 1")
    'A.WriteLine ("   LOAD = 1")
    A.WriteLine ("   VECTOR(SORT1,REAL)=ALL")
    A.WriteLine ("   SPCFORCES(SORT1,REAL)=ALL")
    A.WriteLine ("   MEFFMASS(SUMMARY,MEFFM,FRACSUM,PARTFAC)")
    A.WriteLine ("")
    A.WriteLine ("")
    A.WriteLine ("BEGIN BULK")
    A.WriteLine ("PARAM    POST    -1")
    A.WriteLine ("PARAM    AUTOSPC NO")
    A.WriteLine ("PARAM    K6ROT  10.")
    A.WriteLine ("PARAM    WTMASS .001")
    A.WriteLine ("PARAM    SNORM  30.")
    A.WriteLine ("PARAM   PRTMAXIM YES")
    A.WriteLine ("EIGRL    1              2000.    10      0")
    A.WriteLine ("")
    A.WriteLine ("")
    
    'Nodes
    For i = 0 To Altura
        A.WriteLine (Integer_to_Nastran("GRID") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran("") & Real_to_Nastran("0") & Real_to_Nastran(CStr(i)) & Real_to_Nastran("0"))
    Next i
    
    If Masa1 = "SI" Then
        A.WriteLine ("$Nodo para la Masa 1")
        A.WriteLine (Integer_to_Nastran("GRID") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran("") & Real_to_Nastran(CStr(0 + X1)) & Real_to_Nastran(CStr(i - 1 + Y1)) & Real_to_Nastran("0"))
    End If
    If Masa2 = "SI" Then
        A.WriteLine ("$Nodo para la Masa 2")
        A.WriteLine (Integer_to_Nastran("GRID") & Integer_to_Nastran(CStr(i + 2)) & Integer_to_Nastran("") & Real_to_Nastran(CStr(0 + X2)) & Real_to_Nastran(CStr(i - 1 + Y2)) & Real_to_Nastran("0"))
    End If
    
    
    'Materials
    A.WriteLine ("")
    A.WriteLine ("")
    A.WriteLine ("$Real bracket material")
    A.WriteLine (Integer_to_Nastran("MAT1") & Integer_to_Nastran("1") & Real_to_Nastran(CStr(Young)) & Real_to_Nastran("") & Real_to_Nastran(CStr(Poisson)) & Real_to_Nastran(CStr(Densidad)))
    A.WriteLine ("$Material for the beam-link of the punctual mass")
    A.WriteLine (Integer_to_Nastran("MAT1") & Integer_to_Nastran("2") & Real_to_Nastran(CStr(7000000)) & Real_to_Nastran("") & Real_to_Nastran(CStr(Poisson)) & Real_to_Nastran(CStr(0.0000000001)))
    
    'Properties
    A.WriteLine ("")
    A.WriteLine ("")
    A.WriteLine ("$Real bracket property")
    A.WriteLine (Integer_to_Nastran("PBARL") & Integer_to_Nastran("1") & Integer_to_Nastran("1") & Integer_to_Nastran("") & Integer_to_Nastran("ROD"))
    A.WriteLine (Integer_to_Nastran("") & Real_to_Nastran(CStr(Diametro / 2)) & Real_to_Nastran("0"))
    A.WriteLine ("$Property for the beam-link of the punctual mass")
    A.WriteLine (Integer_to_Nastran("PBARL") & Integer_to_Nastran("2") & Integer_to_Nastran("2") & Integer_to_Nastran("") & Integer_to_Nastran("ROD"))
    A.WriteLine (Integer_to_Nastran("") & Real_to_Nastran(CStr(Diametro * 10)) & Real_to_Nastran("0"))
    
    
    'CBAR Elements
    A.WriteLine ("")
    A.WriteLine ("")
    For i = 0 To Altura - 1
        A.WriteLine (Integer_to_Nastran("CBAR") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran("1") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran(CStr(i + 2)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(1)))
    Next i
    
    If Masa1 = "SI" Then
        A.WriteLine ("$Barra para la Masa 1")
        A.WriteLine (Integer_to_Nastran("CBAR") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran("2") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran(CStr(i + 2)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(1)))
    End If
    If Masa2 = "SI" Then
        A.WriteLine ("$Barra para la Masa 2")
        A.WriteLine (Integer_to_Nastran("CBAR") & Integer_to_Nastran(CStr(i + 2)) & Integer_to_Nastran("2") & Integer_to_Nastran(CStr(i + 1)) & Integer_to_Nastran(CStr(i + 3)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(0)) & Real_to_Nastran(CStr(1)))
    End If
    
    If Masa1 = "SI" Then
        A.WriteLine ("$Masa puntual 1")
        A.WriteLine (Integer_to_Nastran("CONM2") & Integer_to_Nastran(CStr(i + 3)) & Integer_to_Nastran(CStr(i + 2)) & Integer_to_Nastran("0") & Real_to_Nastran(CStr(M1)) & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0"))
        A.WriteLine (Real_to_Nastran("") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0"))
    End If
    If Masa2 = "SI" Then
        A.WriteLine ("$Masa puntual 2")
        A.WriteLine (Integer_to_Nastran("CONM2") & Integer_to_Nastran(CStr(i + 4)) & Integer_to_Nastran(CStr(i + 3)) & Integer_to_Nastran("0") & Real_to_Nastran(CStr(M2)) & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0"))
        A.WriteLine (Real_to_Nastran("") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0") & Real_to_Nastran("0"))
    End If
    
    
    'Boundary conditions
    A.WriteLine ("")
    A.WriteLine ("")
    A.WriteLine (Integer_to_Nastran("SPC1") & Integer_to_Nastran("1") & Integer_to_Nastran("123456") & Integer_to_Nastran("1"))
    
    
    A.WriteLine ("")
    A.WriteLine ("")
    A.WriteLine ("ENDDATA")
    A.Close
    '-----------------------------------------------------
    
End Sub
