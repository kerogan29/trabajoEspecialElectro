''' Trabajo especial de Electrómagnetismo
    Resolución mediante el método de de los momentos
    Facultad de Ingeniería - Universidad de Montevideo
    2016
'''
# Modules goes here
import numpy as np
from scipy import constants as s

# Code starts here
class Tubo:                 # Clase que representa a la seccion transversal del tubo
                            # Definir varibles
    def __init__(self, n):  # Constructor; n -> numero de bloques por lado del tubo cuadrado
        self.gridSize=n     # gridsize -> Lado del tubo
        if n % 2 == 0:      # solo puede ser par (tenemos a/2)
            self.seccTubo = np.zeros((n, n)) # Definir seccion tranversal
        else:
            raise ErrorNum("El numero de divisiones debe ser par") # Levanta una excepción si el numero no es par

    def setVoltage(self,nrocond,value): # Se indica cual conductor estara a Q=1c
        n=self.gridSize
        if nrocond==1:      # Cond 1
            for c in range(0,int(n/2)):
                self.seccTubo[0][c]=value
        elif nrocond==2:    # Cond 2
            for c in range(0, int(n/2)-1):
                self.seccTubo[0][int(n/2)+c]= value
            for c in range(0, n):
                self.seccTubo[c][n-1]= value
        elif nrocond==3:    # Cond 3
             for c in range(0,int(n/2)-1):
                 self.seccTubo[int(n/2)+c][0]= value
             for c in range(0,n):
                 self.seccTubo[n-1][c]= value
        elif nrocond==4:    # Cond 4
            for c in range(0,int(n/2)):
                self.seccTubo[c][0]= value

    def matrixDistAux(self):
        ret=np.zeros((self.gridSize,self.gridSize), dtype=(np.int,(2,2)))
        i=0
        while i<ret.shape[0]:
            j=0
            while j<ret.shape[0]:
                ret[i][j]= (i,j) # Tupla de valores
                j+=1
            i+=1
        return ret

    def lMatrix(self):
        n=self.gridSize
        matdist= self.matrixDistAux()
        l=np.zeros((2*n,2*n)) # Defino la matriz l como una cuadrada de n/2+n
        ratio=0.282 # valor usado cuando i=j
        i=0 # Fila
        while i< l.shape[0]:
            j=0 # Columnna
            while j<l.shape[0]:
                if i==j :
                    temp=ratio*np.log(1+np.sqrt(2))/(2*s.pi)
                    l[i][j]= temp # FALTA EPSILON_0
                elif (i<n/2 and j<n/2):
                    l[i][j]=1/(s.pi*np.abs(i-j)) # FALTA EPSILON_0
                    l[j][i]=l[i][j]
                elif (i>n/2 and j>n/2):
                    if ((i-n/2)<n and (j-n/2)<n): # Me encuntro en segmento horizontal del conductor 3
                        l[i][j] = 1 / (s.pi * np.abs(i - j))  # FALTA EPSILON_0
                        l[j][i] = l[i][j]
                    elif ((i-n/2)>n and (j-n/2)>n):
                        l[i][j]= 1 / (s.pi * np.abs(i - j))
                        l[j][i]= l[i][j]
                elif (i<n/2 and j>n/2):
                    if j<(3/2)*n:
                        l[i][j]=1/(s.pi * np.sqrt(np.power((matdist[0][i][0][0]-matdist[n-1][j-n/2][0][0]),2)+np.power((matdist[0][i][0][1]-matdist[n-1][j-n/2][0][1]),2)))
                        l[j][i]= l[i][j]
                    else:
                        l[i][j] = 1/(s.pi * np.sqrt(np.power((matdist[0][i][0][0] - matdist[j - (3/2)*n][0][0][0]), 2)+ np.power((matdist[0][i][0][1] - matdist[j - (3/2)*n][0][0][1]), 2)))
                        l[j][i] = l[i][j]
                j+=1
            i+=1
        return l

class ErrorNum(BaseException):
    def __init__(self, mensaje):
        print(mensaje)

if __name__ == "__main__":
    t=Tubo(4)
    print(t.lMatrix())
    #print("Mesh Relaxation method")
    #print("Info:\n\tTamaño de Mesh:",t.gridSize)
    #cij=t.capacitanceCoeff(1,3,0.00001)
    #print("\tC13=\t",cij)
    #l=Tubo(10)
    #cji=l.capacitanceCoeff(3,1,0.00001)
    #print("\tC31=\t",cji)
    #print("\tPromedio(C13+C31/2)= ",(cij+cji)/2)
    #print("\tdiff=\t",np.abs(cij-cji))
