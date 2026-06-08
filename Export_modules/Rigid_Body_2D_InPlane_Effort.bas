Attribute VB_Name = "Rigid_Body_2D_InPlane_Effort"
'------------------------------------------------------------------------------
'
' Module de Calcul de Distribution d'un Corps Rigide sur Fondation Elastique
' 2D Dans Plan (Plan XY: FX, FY et MZ)
'
' Version 1.0 ( 26/11/04 )
'
'------------------------------------------------------------------------------
' Liste des Fonctions Définies dans le Module
'------------------------------------------------------------------------------

' - Création Matrice de Rigidité 2D d'un Système Corps Rigide sur Ressorts: Dans le Plan
' Rigid_Body_2D_InPlane_Stiffness_Matrix(Efficiency, Point_Location, Stiff_Direction_1_AnglevsXaxis_Deg,
' Stiff_Weighting, NB_Points)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points*2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points), Stiff_Weighting(NB_Points * 2 Directions)
' Sortie Matricielle: Matrice Rigidité à l'Origine du Repère (3*3)

' - Distribution d'Efforts d'un Système Corps Rigide 2D sur Ressorts: Dans le Plan
' Rigid_Body_2D_InPlane_Effort_Distribution(Efficiency, Point_Location, Stiff_Direction_1_AnglevsXaxis_Deg,
' Stiff_Weighting, NB_Points, Loading_Point, Loading_Effort)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points*2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points), Stiff_Weighting(NB_Points)
' Sortie Matricielle: Matrice des Efforts Distribués Dans le Plan (NB_Points * 2 Coordonnées) dans le
' repère local des points

' - Sommation d'Efforts Locaux défini dans un Repère Global
' Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord(Efficiency, Point_Location,
'  Local_Loading, NB_Points, Load_Resultant_Point)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points * 2 Coordonnées)
' Local_Loading(NB_Points * 1 Coordonnées), Load_Resultant_Point(2 Coordonnées)
' Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)

' - Sommation d'Efforts Locaux définis dans un Repère Local
' Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord(Efficiency, Point_Location,
'  Stiff_Direction_1_AnglevsXaxis_Deg, Local_Loading, NB_Points, Load_Resultant_Point)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points * 2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points * 1 Coordonnée)
' Local_Loading(NB_Points * 2 Coordonnées), Load_Resultant_Point(2 Coordonnées)
' Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)

' - Changement de Base d'Expression d'Efforts
' Rigid_Body_2D_InPlane_Effort_Coord_Change(Direction_1_AnglevsXaxis_Deg, Set_Of_Effort)
' Sortie Matricielle: Matrice des Efforts dans Nouveau Repère (1 Ligne * 2 Composantes)

'------------------------------------------------------------------------------
Option Explicit
'------------------------------------------------------------------------------

'******************************************************************************
' Initialisation: Informations sur Contenu des Fonctions
'******************************************************************************
Sub Rigid_Body_2D_InPlane_Effort_Initialize()
  
' Rigid_Body_2D_InPlane_Stiffness_Matrix
Application.MacroOptions Macro:="Rigid_Body_2D_InPlane_Stiffness_Matrix", Description:= _
    "Création Matrice de Rigidité 2D d'un Système Corps Rigide sur Ressorts: Dans le Plan" & vbCrLf _
    & "Entrées Matricielles: Tableaux(NB_Points * X Coord)" & vbCrLf _
    & "X = 1 ou 2 selon donnée" & vbCrLf _
    & "Sortie Matricielle: Matrice Rigidité à l'Origine du Repère (3*3)"

' Rigid_Body_2D_InPlane_Effort_Distribution
Application.MacroOptions Macro:="Rigid_Body_2D_InPlane_Effort_Distribution", Description:= _
    "Distribution d'Efforts d'un Système Corps Rigide 2D sur Ressorts: Dans le Plan" & vbCrLf _
    & "Entrées Matricielles: Tableau(NB_Points * X Coordonnées)" _
    & "(X = 1 ou 2 selon donnée)" & vbCrLf _
    & "Sortie Matricielle: Matrice des Efforts Distribués (NB_Points * 2 Coordonnées) Repère Local"

' Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord
Application.MacroOptions Macro:="Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord", Description:= _
    "Sommation d'Efforts Locaux Dans le Plan Définis dans un Repère Global" & vbCrLf _
    & "Entrées Matricielles: Tableaux(NB_Points * X Coord)" _
    & " Point_Location, Local_Loading: X=2" & vbCrLf _
    & "Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)"

' Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord
Application.MacroOptions Macro:="Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord", Description:= _
    "Sommation d'Efforts Locaux Dans le Plan  définis dans un Repère Local" & vbCrLf _
    & "Entrées Matricielles: Tableaux(NB_Points * X Coord)" _
    & " Point_Location, Local_Loading: X=2" & vbCrLf _
    & "Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)"

' Rigid_Body_2D_InPlane_Effort_Coord_Change
Application.MacroOptions Macro:="Rigid_Body_2D_InPlane_Effort_Coord_Change", Description:= _
    "Changement de Base d'Expression d'Efforts" & vbCrLf _
    & "Sortie Matricielle: Matrice des Efforts dans Nouveau Repère (1 Ligne * 2 Composantes)"

End Sub


'******************************************************************************
' Création Matrice de Rigidité 2D d'un Système Corps Rigide sur Ressorts: Dans le Plan
' Rigid_Body_2D_InPlane_Stiffness_Matrix(Efficiency, Point_Location, Stiff_Direction_1_AnglevsXaxis_Deg,
' Stiff_Weighting, NB_Points)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points*2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points), Stiff_Weighting(NB_Points * 2 Directions)
' Sortie Matricielle: Matrice Rigidité à l'Origine du Repère (3*3)
'******************************************************************************
Function Rigid_Body_2D_InPlane_Stiffness_Matrix(Efficiency As Variant, _
  Point_Location As Variant, Stiff_Direction_1_AnglevsXaxis_Deg As Variant, Stiff_Weighting As Variant, _
  NB_Points As Integer) As Variant
Attribute Rigid_Body_2D_InPlane_Stiffness_Matrix.VB_Description = "Création Matrice de Rigidité 2D d'un Système Corps Rigide sur Ressorts: Dans le Plan\r\nEntrées Matricielles: Tableaux(NB_Points * X Coord)\r\nX = 1 ou 2 selon donnée\r\nSortie Matricielle: Matrice Rigidité à l'Origine du Repère (3*3)"
Attribute Rigid_Body_2D_InPlane_Stiffness_Matrix.VB_ProcData.VB_Invoke_Func = " \n14"

Dim Ind_Disp As Integer, Ind_Point As Integer, Ind_DDL As Integer
Dim tmp_Global_Disp(1 To 3, 1 To 3) As Double, tmp_Unitary_Global_Disp(1 To 3, 1 To 1) As Double
Dim tmp_Local_Disp_RG(1 To 2, 1 To 1) As Double
Dim tmp_Local_Disp_RL(1 To 2, 1 To 1) As Double
Dim tmp_Local_Effort_RG(1 To 2, 1 To 1) As Double
Dim tmp_Local_Effort_RL(1 To 2, 1 To 1) As Double
Dim tmp_Stiffness_Matrix(1 To 3, 1 To 3) As Double
Dim tmp_Angle_Rad As Double

' Constitution de la Matrice de Rigidité du Système: Vecteurs Déplacement Unitaire
For Ind_Disp = 1 To 3
  tmp_Global_Disp(Ind_Disp, Ind_Disp) = 1
Next Ind_Disp

' Parcours par Point de Rigidité
For Ind_Point = 1 To NB_Points
  If Efficiency(Ind_Point) = 1 Then
' Parcours par Déplacement Unitaire
    For Ind_Disp = 1 To 3
' Vecteur Déplacement Unitaire
      For Ind_DDL = 1 To 3
        tmp_Unitary_Global_Disp(Ind_DDL, 1) = tmp_Global_Disp(Ind_DDL, Ind_Disp)
      Next Ind_DDL
' Vecteur Déplacement Local au Point: Expression dans Repère Global
      tmp_Local_Disp_RG(1, 1) = tmp_Unitary_Global_Disp(1, 1) _
        - tmp_Unitary_Global_Disp(3, 1) * Point_Location(Ind_Point, 2)
      tmp_Local_Disp_RG(2, 1) = tmp_Unitary_Global_Disp(2, 1) _
        + tmp_Unitary_Global_Disp(3, 1) * Point_Location(Ind_Point, 1)
' Vecteur Déplacement Local au Point: Expression dans Repère Local
      tmp_Angle_Rad = Stiff_Direction_1_AnglevsXaxis_Deg(Ind_Point, 1) / 180 _
        * Application.pi()
      tmp_Local_Disp_RL(1, 1) = tmp_Local_Disp_RG(1, 1) * Cos(tmp_Angle_Rad) _
        + tmp_Local_Disp_RG(2, 1) * Sin(tmp_Angle_Rad)
      tmp_Local_Disp_RL(2, 1) = -tmp_Local_Disp_RG(1, 1) * Sin(tmp_Angle_Rad) _
        + tmp_Local_Disp_RG(2, 1) * Cos(tmp_Angle_Rad)
' Effort Local Ressort: Expression dans Repère Local
      tmp_Local_Effort_RL(1, 1) = tmp_Local_Disp_RL(1, 1) * Stiff_Weighting(Ind_Point, 1)
      tmp_Local_Effort_RL(2, 1) = tmp_Local_Disp_RL(2, 1) * Stiff_Weighting(Ind_Point, 2)
' Effort Local Ressort: Expression dans Repère Global
      tmp_Local_Effort_RG(1, 1) = tmp_Local_Effort_RL(1, 1) * Cos(tmp_Angle_Rad) _
        - tmp_Local_Effort_RL(2, 1) * Sin(tmp_Angle_Rad)
      tmp_Local_Effort_RG(2, 1) = tmp_Local_Effort_RL(1, 1) * Sin(tmp_Angle_Rad) _
        + tmp_Local_Effort_RL(2, 1) * Cos(tmp_Angle_Rad)
' Effort Ressort: Expression dans Repère Global à l'Origine du Repère
' Constitution Matrice de Rigidité du Système
      tmp_Stiffness_Matrix(1, Ind_Disp) = tmp_Stiffness_Matrix(1, Ind_Disp) _
        + tmp_Local_Effort_RG(1, 1)
      tmp_Stiffness_Matrix(2, Ind_Disp) = tmp_Stiffness_Matrix(2, Ind_Disp) _
        + tmp_Local_Effort_RG(2, 1)
      tmp_Stiffness_Matrix(3, Ind_Disp) = tmp_Stiffness_Matrix(3, Ind_Disp) _
        - tmp_Local_Effort_RG(1, 1) * Point_Location(Ind_Point, 2) _
        + tmp_Local_Effort_RG(2, 1) * Point_Location(Ind_Point, 1)
    Next Ind_Disp
  End If
Next Ind_Point

' Sortie Matrice de Rigidité du Système
Rigid_Body_2D_InPlane_Stiffness_Matrix = tmp_Stiffness_Matrix

End Function


'******************************************************************************
' Distribution d'Efforts d'un Système Corps Rigide 2D sur Ressorts: Dans le Plan
' Rigid_Body_2D_InPlane_Effort_Distribution(Efficiency, Point_Location, Stiff_Direction_1_AnglevsXaxis_Deg,
' Stiff_Weighting, NB_Points, Loading_Point, Loading_Effort)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points*2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points), Stiff_Weighting(NB_Points)
' Sortie Matricielle: Matrice des Efforts Distribués Dans le Plan (NB_Points * 2 Coordonnées) dans le
' repère local des points
'******************************************************************************
Function Rigid_Body_2D_InPlane_Effort_Distribution(Efficiency As Variant, _
  Point_Location As Variant, Stiff_Direction_1_AnglevsXaxis_Deg As Variant, _
  Stiff_Weighting As Variant, NB_Points As Integer, _
  Loading_Point As Variant, Loading_Effort As Variant) As Variant
Attribute Rigid_Body_2D_InPlane_Effort_Distribution.VB_Description = "Distribution d'Efforts d'un Système Corps Rigide 2D sur Ressorts: Dans le Plan\r\nEntrées Matricielles: Tableau(NB_Points * X Coordonnées)(X = 1 ou 2 selon donnée)\r\nSortie Matricielle: Matrice des Efforts Distribués (NB_Points * 2 Coordonnées) Repère Local"
Attribute Rigid_Body_2D_InPlane_Effort_Distribution.VB_ProcData.VB_Invoke_Func = " \n14"

Dim tmp_Stiffness_Matrix As Variant, tmp_Inv_Stiffness_Matrix As Variant
Dim tmp_Stiffness_Determinant As Double
Dim tmp_Effort_Distribue() As Double
Dim tmp_Deformation_Globale_Var As Variant
Dim tmp_Torseur_Origine(1 To 3, 1 To 1) As Double
Dim Ind_Point As Integer
Dim tmp_Local_Disp_RG(1 To 2, 1 To 1) As Double
Dim tmp_Local_Disp_RL(1 To 2, 1 To 1) As Double
Dim tmp_Local_Effort_RL(1 To 2, 1 To 1) As Double
Dim tmp_Angle_Rad As Double

ReDim tmp_Effort_Distribue(1 To NB_Points, 1 To 2)

' Obtention Matrice de Rigidité du Système
tmp_Stiffness_Matrix = Rigid_Body_2D_InPlane_Stiffness_Matrix(Efficiency, Point_Location, _
  Stiff_Direction_1_AnglevsXaxis_Deg, Stiff_Weighting, NB_Points)

' Calcul du Déterminant de la Matrice de Rigidité
tmp_Stiffness_Determinant = Application.MDeterm(tmp_Stiffness_Matrix)
' Test Inversibilité Matrice
If tmp_Stiffness_Determinant = 0 Then
  For Ind_Point = 1 To NB_Points
    tmp_Effort_Distribue(Ind_Point, 1) = "Singular Matrix"
    tmp_Effort_Distribue(Ind_Point, 2) = "Singular Matrix"
  Next Ind_Point
  Rigid_Body_2D_InPlane_Effort_Distribution = tmp_Effort_Distribue
  Exit Function
End If

' Inversion Matrice de Rigidité
tmp_Inv_Stiffness_Matrix = Application.MInverse(tmp_Stiffness_Matrix)
' Expression du Torseur des Efforts Appliqués à l'Origine du Repère
tmp_Torseur_Origine(1, 1) = Loading_Effort(1, 1)
tmp_Torseur_Origine(2, 1) = Loading_Effort(2, 1)
tmp_Torseur_Origine(3, 1) = Loading_Effort(3, 1) _
  - Loading_Effort(1, 1) * Loading_Point(2, 1) + Loading_Effort(2, 1) * Loading_Point(1, 1)
' Détermination Déformation dans le Repère Global
tmp_Deformation_Globale_Var = Application.MMult(tmp_Inv_Stiffness_Matrix, tmp_Torseur_Origine)

' Parcours par Point de Rigidité
For Ind_Point = 1 To NB_Points
  If Efficiency(Ind_Point) = 1 Then
' Détermination Déformation Locale: Repère Global
    tmp_Local_Disp_RG(1, 1) = tmp_Deformation_Globale_Var(1, 1) _
      - tmp_Deformation_Globale_Var(3, 1) * Point_Location(Ind_Point, 2)
    tmp_Local_Disp_RG(2, 1) = tmp_Deformation_Globale_Var(2, 1) _
      + tmp_Deformation_Globale_Var(3, 1) * Point_Location(Ind_Point, 1)
' Détermination Déformation Locale: Repère Local
    tmp_Angle_Rad = Stiff_Direction_1_AnglevsXaxis_Deg(Ind_Point, 1) / 180 _
      * Application.pi()
    tmp_Local_Disp_RL(1, 1) = tmp_Local_Disp_RG(1, 1) * Cos(tmp_Angle_Rad) _
      + tmp_Local_Disp_RG(2, 1) * Sin(tmp_Angle_Rad)
    tmp_Local_Disp_RL(2, 1) = -tmp_Local_Disp_RG(1, 1) * Sin(tmp_Angle_Rad) _
      + tmp_Local_Disp_RG(2, 1) * Cos(tmp_Angle_Rad)
' Détermination Effort Local: Repère Global
    tmp_Effort_Distribue(Ind_Point, 1) = tmp_Local_Disp_RL(1, 1) * Stiff_Weighting(Ind_Point, 1)
    tmp_Effort_Distribue(Ind_Point, 2) = tmp_Local_Disp_RL(2, 1) * Stiff_Weighting(Ind_Point, 2)
  End If
Next Ind_Point

Rigid_Body_2D_InPlane_Effort_Distribution = tmp_Effort_Distribue

End Function


'******************************************************************************
' Sommation d'Efforts Locaux défini dans un Repère Global
' Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord(Efficiency, Point_Location,
'  Local_Loading, NB_Points, Load_Resultant_Point)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points * 2 Coordonnées)
' Local_Loading(NB_Points * 2 Coordonnées), Load_Resultant_Point(2 Coordonnées)
' Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)
'******************************************************************************
Function Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord(Efficiency As Variant, _
  Point_Location As Variant, Local_Loading As Variant, NB_Points As Integer, _
  Load_Resultant_Point As Variant) As Variant
Attribute Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord.VB_Description = "Sommation d'Efforts Locaux Dans le Plan Définis dans un Repère Global\r\nEntrées Matricielles: Tableaux(NB_Points * X Coord) Point_Location, Local_Loading: X=2\r\nSortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)"
Attribute Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord.VB_ProcData.VB_Invoke_Func = " \n14"

Dim tmp_Torseur_Resultant(1 To 3, 1 To 1) As Double
Dim Ind_Point As Integer, Ind_Coord As Integer

' Parcours par Point de Rigidité
For Ind_Point = 1 To NB_Points
  If Efficiency(Ind_Point) = 1 Then
    tmp_Torseur_Resultant(1, 1) = tmp_Torseur_Resultant(1, 1) + Local_Loading(Ind_Point, 1)
    tmp_Torseur_Resultant(2, 1) = tmp_Torseur_Resultant(2, 1) + Local_Loading(Ind_Point, 2)
    tmp_Torseur_Resultant(3, 1) = tmp_Torseur_Resultant(3, 1) _
      - Local_Loading(Ind_Point, 1) * (Point_Location(Ind_Point, 2) - Load_Resultant_Point(2)) _
      + Local_Loading(Ind_Point, 2) * (Point_Location(Ind_Point, 1) - Load_Resultant_Point(1)) _
      + Local_Loading(Ind_Point, 3)
  End If
Next Ind_Point

Rigid_Body_2D_InPlane_Effort_Resultant_GlbCoord = tmp_Torseur_Resultant

End Function


'******************************************************************************
' Sommation d'Efforts Locaux définis dans un Repère Local
' Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord(Efficiency, Point_Location,
'  Stiff_Direction_1_AnglevsXaxis_Deg, Local_Loading, NB_Points, Load_Resultant_Point)
' Entrée Matricielle: Efficiency(NB_Points), Point_Location(NB_Points * 2 Coordonnées)
' Stiff_Direction_1_AnglevsXaxis_Deg(NB_Points * 1 Coordonnée)
' Local_Loading(NB_Points * 2 Coordonnées), Load_Resultant_Point(2 Coordonnées)
' Sortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)
'******************************************************************************
Function Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord(Efficiency As Variant, _
  Point_Location As Variant, Stiff_Direction_1_AnglevsXaxis_Deg As Variant, _
  Local_Loading As Variant, NB_Points As Integer, _
  Load_Resultant_Point As Variant) As Variant
Attribute Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord.VB_Description = "Sommation d'Efforts Locaux Dans le Plan  définis dans un Repère Local\r\nEntrées Matricielles: Tableaux(NB_Points * X Coord) Point_Location, Local_Loading: X=2\r\nSortie Matricielle: Matrice du Torseur Résultant (3 Composantes * 1 Colonne)"
Attribute Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord.VB_ProcData.VB_Invoke_Func = " \n14"

Dim tmp_Torseur_Resultant(1 To 3, 1 To 1) As Double
Dim tmp_Torseur_RG As Variant
Dim tmp_Torseur_RL(1 To 2) As Double
Dim Ind_Point As Integer, Ind_Coord As Integer

' Parcours par Point de Rigidité
For Ind_Point = 1 To NB_Points
  If Efficiency(Ind_Point) = 1 Then
    tmp_Torseur_RL(1) = Local_Loading(Ind_Point, 1)
    tmp_Torseur_RL(2) = Local_Loading(Ind_Point, 2)
    tmp_Torseur_RG = Rigid_Body_2D_InPlane_Effort_Coord_Change( _
      -Stiff_Direction_1_AnglevsXaxis_Deg(Ind_Point, 1), tmp_Torseur_RL)
    tmp_Torseur_Resultant(1, 1) = tmp_Torseur_Resultant(1, 1) + tmp_Torseur_RG(1)
    tmp_Torseur_Resultant(2, 1) = tmp_Torseur_Resultant(2, 1) + tmp_Torseur_RG(2)
    tmp_Torseur_Resultant(3, 1) = tmp_Torseur_Resultant(3, 1) _
      - tmp_Torseur_RG(1) * (Point_Location(Ind_Point, 2) - Load_Resultant_Point(2)) _
      + tmp_Torseur_RG(2) * (Point_Location(Ind_Point, 1) - Load_Resultant_Point(1))
  End If
Next Ind_Point

Rigid_Body_2D_InPlane_Effort_Resultant_LclCoord = tmp_Torseur_Resultant

End Function



'******************************************************************************
' Changement de Base d'Expression d'Efforts
' Rigid_Body_2D_InPlane_Effort_Coord_Change(Direction_1_AnglevsXaxis_Deg, Set_Of_Effort)
' Sortie Matricielle: Matrice des Efforts dans Nouveau Repère (1 Ligne * 2 Composantes)
'******************************************************************************
Function Rigid_Body_2D_InPlane_Effort_Coord_Change(Direction_1_AnglevsXaxis_Deg As Double, _
  Set_Of_Effort As Variant) As Variant
Attribute Rigid_Body_2D_InPlane_Effort_Coord_Change.VB_Description = "Changement de Base d'Expression d'Efforts\r\nSortie Matricielle: Matrice des Efforts dans Nouveau Repère (1 Ligne * 2 Composantes)"
Attribute Rigid_Body_2D_InPlane_Effort_Coord_Change.VB_ProcData.VB_Invoke_Func = " \n14"
  
Dim tmp_Angle_Rad As Double
Dim tmp_Torseur_New_Coord(1 To 2) As Double

tmp_Angle_Rad = Direction_1_AnglevsXaxis_Deg / 180 * Application.pi()
tmp_Torseur_New_Coord(1) = Set_Of_Effort(1) * Cos(tmp_Angle_Rad) + Set_Of_Effort(2) * Sin(tmp_Angle_Rad)
tmp_Torseur_New_Coord(2) = -Set_Of_Effort(1) * Sin(tmp_Angle_Rad) + Set_Of_Effort(2) * Cos(tmp_Angle_Rad)

Rigid_Body_2D_InPlane_Effort_Coord_Change = tmp_Torseur_New_Coord

End Function


