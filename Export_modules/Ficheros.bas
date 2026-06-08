Attribute VB_Name = "Ficheros"
Function nombre_libro()
    ' Devuelve el nombre del libro actual
    Dim nombre As String
    nombre = ThisWorkbook.Name
    nombre_libro = nombre
End Function

Function ruta_libro()
    ' Devuelve la ruta del libro actual
    Dim nombre As String
    nombre = ThisWorkbook.path
    ruta_libro = nombre
End Function

Function peso_archivo(path As String, nombre As String)
    ' Devuelve el tamaño de un archivo en Bytes
    Ruta = path & "\" & nombre
    peso_archivo = FileLen(Ruta)
End Function

Function crea_carpeta(path As String, nombre As String)
    ' Crea una carpeta
    ' Abajo hay hotra hecha con objetos
    Ruta = path & "\" & nombre
    MkDir (Ruta)
End Function

Sub cmd()
    ' Abre el símbolo del sistema (consola de MS-DOS)
    dblRetVal = Shell("cmd.exe")
End Sub

Sub Importar_txt()
    Dim hoja_nombre As String
    Dim libro_nombre As String
    Dim libro_auxiliar_nombre As String
    hoja_nombre = ActiveSheet.Name
    libro_nombre = nombre_libro()

    ' Se abre el ".txt" en un libro auxiliar
    ' En la siguiente línea hay que dar el nombre y ruta del archivo de texto que se quiere importar
    Workbooks.OpenText Filename:="C:\Users\AndrésJavier\Desktop\KK\kkk.txt", _
        Origin:=xlMSDOS, StartRow:=1, DataType:=xlDelimited, TextQualifier:= _
        xlDoubleQuote, ConsecutiveDelimiter:=False, Tab:=False, Semicolon:=False _
        , Comma:=False, Space:=False, Other:=False, FieldInfo:=Array(1, 1), _
        TrailingMinusNumbers:=True
    ' Se almacena el nombre del libro auxiliar
    libro_auxiliar_nombre = ActiveWorkbook.Name


    ' Se selecciona y copia toda la columna
    Columns(1).Select
    Selection.Copy
    
    ' Se activa el libro original
    Windows(libro_nombre).Activate
    Range("A1").Select
    
    ' Se pega en la hoja desde la que se llamó a la macro
    ActiveSheet.Paste
    Application.CutCopyMode = False
    
    'Se cierra el libro auxiliar
    Windows(libro_auxiliar_nombre).Activate
    ActiveWindow.Close
    Range("A1").Select

End Sub
Sub OpenWeb()
    ' Abre el Internet Explorer y lo direcciona a una web
    ' CreateObject Method creates an Automation object of the specified class
    Set Web = CreateObject("InternetExplorer.Application")
    Web.Visible = True
    Web.Navigate "www.microsoft.com"
End Sub
Sub crea_carpeta2(path As String)
    ' Por ejemplo, para crear una carpeta llamando a esta subrutina: call crea_carpeta2("D:\KK")
    On Error GoTo error:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FOLDER = FSO.CreateFolder(path)
error:
End Sub
Function existe_carpeta(path As String)
    ' Devuelve True o False si la carpeta existe o no
    Set FSO = CreateObject("Scripting.FileSystemObject")
    existe_carpeta = FSO.FolderExists(path)
End Function
Sub borra_carpeta(path As String)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FOLDER_TO_DELETE = FSO.GetFolder(path)
    FOLDER_TO_DELETE.Delete
End Sub

Sub txt()
    ' Crea un archivo de texto y lo edita
    ' The CreateObject function returns the FileSystemObject (fso)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    ' The CreateTextFile method then creates the file as a TextStream object (a)
    Set A = FSO.createtextfile("D:\kk.txt", True) 'object.CreateTextFile(filename[, overwrite[, unicode]]) --> En este caso, "objet" es el objeto "fso" creado anteriormente
    ' And the WriteLine method writes a line of text to the created text file
    A.WriteLine ("This is a test.")
    A.WriteLine ("This is a test2.")
    ' The Close method flushes the buffer and closes the file
    A.Close
End Sub

Function existe_archivo(path As String)
    ' Devuelve True o False si el archivo existe o no
    Set FSO = CreateObject("Scripting.FileSystemObject")
    existe_archivo = FSO.FileExists(path)
End Function
Sub borra_archivo(path As String)
    ' Borra el archivo especificado
    ' Por ejemplo, para borrar un archivo llamando a esta subrutina: call borra_archivo("D:\kk.txt")
    On Error GoTo error:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FILE_TO_DELETE = FSO.GetFile(path)
    FILE_TO_DELETE.Delete
    ' Todo el proceso se podría hacer en una sola línea sustituyendo:     Set FILE_TO_DELETE = CreateObject("Scripting.FileSystemObject").getfile(path).Delete
error:
End Sub





'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''' All the functions bellow are used at some point along the program. Their names are quite explanatory'''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'Estas funciones no son mías. Están imitadas más arriba

Function FOLDER_SIZE(ByVal FOLDER_PATH As String) As Variant
    On Error GoTo err:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FOLDER = FSO.GetFolder(FOLDER_PATH)
    FOLDER_SIZE = FOLDER.Size
err:
End Function

Function FILE_SIZE(ByVal FILE_PATH As String) As Variant
    On Error GoTo err:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set file = FSO.GetFile(FILE_PATH)
    FILE_SIZE = file.Size
err:
End Function

Function FILE_EXIST(file As String)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    FILE_EXIST = FSO.FileExists(file)
End Function

Function FOLDER_EXIST(FOLDER As String)
    Set FSO = CreateObject("Scripting.FileSystemObject")
    FOLDER_EXIST = FSO.FolderExists(FOLDER)
End Function

Sub FILE_DELETE(ByVal file As String)
    On Error GoTo err:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FILE_TO_DELETE = FSO.GetFile(file)
    FILE_TO_DELETE.Delete
err:
End Sub

Sub FOLDER_CREATE(ByVal FOLDER_PATH As String)
    On Error GoTo err:
    Set FSO = CreateObject("Scripting.FileSystemObject")
    Set FOLDER = FSO.CreateFolder(FOLDER_PATH)
err:
End Sub

