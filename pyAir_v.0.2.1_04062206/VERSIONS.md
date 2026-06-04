# ================================ #
#                                  #              
#   VERSION [0.1.0] - 26-03-2026   #
#                                  #  
# ================================ #

### Added
### -----
 - Initial Release


# ================================ #
#                                  #              
#   VERSION [0.2.0] - 26-03-2026   #
#                                  #  
# ================================ #

### Added
### -----
	- ANALYSIS.ARPA
		- El problema de autovalores del sistema de ecuaciones resultante de plantear la 
		resolución del problema por métodos energéticos puede estar resultando en un problema
		mal condicionado. Las señales que pueden estar indicando este problema son:
			- autovalores negativos
			- autovalores muy grandes
			- Modos con oscilaciones absurdas o espureas
		
		  La manera de resolverlo es mediante el metodo de Tikhonov que en este caso consiste en:
			 Q ---> Qlambda = Q +lambda * I
			 T ---> Tlambda = T + lambda * I
