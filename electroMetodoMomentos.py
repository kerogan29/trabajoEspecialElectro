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

    def voltVector(self):
        n=self.gridSize
        vVoltage=np.ones(4*n) # Todos los conductores inicializados en cero
        cte= -2*s.pi*s.epsilon_0
        i=0
        while i< vVoltage.shape[0]:
            if i<= n/2:
                vVoltage[i]=self.seccTubo[0][i]
            elif i>n/2 and i<n:
                vVoltage[i] = self.seccTubo[0][i]
            elif i>=n and i<2*n:
                vVoltage[i] = self.seccTubo[i-n][n-1]
            elif i>=2*n and i<(5/2)*n:
                vVoltage[i] = self.seccTubo[i-int((3/2)*n)][0]
            elif i>=(5/2)*n and i<(7/2)*n:
                vVoltage[i] = self.seccTubo[n-1][i-int((5/2)*n)]
            else:
                vVoltage[i] = self.seccTubo[int(i-(7/2)*n)][n-1]
            i+=1
        return vVoltage*cte

    def matrixDistAux(self):
        ret=[] #np.zeros((self.gridSize,self.gridSize), dtype=(np.int,(2,2)))
        i=0
        while i<self.gridSize:
            j=0
            while j<self.gridSize:
                ret.append((i,j)) # Tupla de valores
                j+=1
            i+=1
        return ret

    def matrixSearch(self,i,j):
        return self.matrixDistAux()[(i*self.gridSize)+j]

    def logDist(self,i1,i2,j1,j2): # Fila 1, Columna 1/ Fila 2, Columna2
        dist1 = np.power(self.matrixSearch(int(i1), int(i2))[0] - self.matrixSearch(int(j1), int(j2))[0], 2)
        dist2 = np.power(self.matrixSearch(int(i1), int(i2))[1] - self.matrixSearch(int(j1), int(j2))[1], 2)
        dist = np.sqrt(dist1 + dist2)
        if dist==0:
            return 0
        else:
            return np.log(dist)

    def llMatrix(self):
        n = self.gridSize
        l = np.zeros((4 * n, 4 * n))  # Defino la matriz l como una cuadrada de n/2+n
        ratio = -(np.log(2)+1)#(1 / 2) * (np.log(1 / 2)) - 1 / 2  # valor usado cuando i=j
        i = 0  # Fila
        while i < l.shape[0]:
            j = i  # Columnna
            while j < l.shape[0]:
                if i == j:              # Interacción del bloque consigo mismo
                    temp = ratio        # ratio*(1/(s.pi*s.epsilon_0))
                    l[i][j] = temp
                elif (i<n/2 and j<n/2): # Lo mismo j<n/2 # Conductor 1 sobre cond1
                    l[i][j] = self.logDist(0,i,0,j)
                    l[j][i] = l[i][j]
                elif (j>=n/2 and j<2*n):
                    if i<n/2:    # Conductor 1 sobre cond2
                        if j<n:
                            l[i][j] = self.logDist(0,i,0,j)
                            l[j][i] = l[i][j]
                        else:
                            l[i][j]= self.logDist(0,i,j-n,n-1)
                            l[j][i]= l[i][j]
                    else:   # Conductor 2 sobre cond2
                        if j<n:   # Segmento de n/2
                            l[i][j] = self.logDist(0, i, 0, j)
                            l[j][i] = l[i][j]
                        else:       #
                            if i<n:
                                l[i][j] = self.logDist(0, i, j-n, n - 1)
                                l[j][i] = l[i][j]
                            else:
                                l[i][j] = self.logDist(i-n,n-1,j-n, n - 1)
                                l[j][i] = l[i][j]
                elif(j>=2*n and j<(7/2)*n): # Cond 3
                    if i<n/2:               # Cond1 sobre 3
                        if j<(5/2)*n:
                            l[i][j] = self.logDist(0,i,j-(2*n),0)
                            l[j][i] = l[i][j]
                        else:
                            l[i][j] = self.logDist(0, i, n-1, j-(5/2)*n)
                            l[j][i] = l[i][j]
                    elif (i>=n/2 and i<2*n): # Dos sobre tres
                        if j<(5/2)*n and i<n:
                            l[i][j] = self.logDist(0, i, j-(2*n),0)
                            l[j][i] = l[i][j]
                        elif j<(5/2)*n and i>=n:
                            l[i][j] = self.logDist(i-n, n-1, j - (2 * n), 0)
                            l[j][i] = l[i][j]
                        elif j>=(5/2)*n and i<n:
                            l[i][j] = self.logDist(0, i, n-1, j-(5/2)*n)
                            l[j][i] = l[i][j]
                        elif j>=(5/2)*n and i>=n:
                            l[i][j] = self.logDist(i-n, n-1, n - 1, j - (5 / 2) * n)
                            l[j][i] = l[i][j]
                    else: # Tres sobre tres
                        if j<(5/2)*n and i<(5/2)*n:
                            l[i][j] = self.logDist(i - (2 * n), 0, j - (2 * n), 0)
                            l[j][i] = l[i][j]
                        elif j>=(5/2)*n and i>=(5/2)*n:
                            l[i][j] = self.logDist(n-1, i-((5/2)*n), n-1, j-((5/2)*n))
                            l[j][i] = l[i][j]
                        else:
                            l[i][j] = self.logDist(i - (2 * n), 0,  n-1, j-(5/2)*n)
                            l[j][i] = l[i][j]
                else: # Cond 4
                    if i<n/2:
                        l[i][j] = self.logDist(0,i, j-(7/2)*n, 0)
                        l[j][i] = l[i][j]
                    elif i>=n/2 and i<n:
                        l[i][j] = self.logDist(0, i, j - (7 / 2) * n, 0)
                        l[j][i] = l[i][j]
                    elif i >= n and i < 2*n:
                        l[i][j] = self.logDist(i-n, n-1, j - (7 / 2) * n, 0)
                        l[j][i] = l[i][j]
                    elif i>=2*n and i<(5/2)*n:
                        l[i][j] = self.logDist(i -2*n, 0, j - (7 / 2) * n, 0)
                        l[j][i] = l[i][j]
                    elif i >= (5 / 2) * n and i < (7 / 2) * n:
                        l[i][j] = self.logDist(n-1,i-(5/2)*n, j - (7 / 2) * n, 0)
                        l[j][i] = l[i][j]
                    else:
                        l[i][j] = self.logDist(i - (7 / 2) * n, 0, j - (7 / 2) * n, 0)
                        l[j][i] = l[i][j]
                j+=1
            i+=1
        return l

    def inv(self,a):
        return np.lianlg(a) # Devuelve la matriz inversa

    def calcCoef(self,i,j):
        n=self.gridSize
        m=self.llMatrix() # Matriz de log de distancias
        if j in (1,2,3,4):
            self.setVoltage(j,1)
        else:
            raise ErrorNum("Coeficiente incorrecto")
        vecCharge=np.dot(m,self.voltVector())
        if i==1:
            c13=0
            for i in range(0,int(n/2)-1):
                c13+=vecCharge[i]
                return np.abs(c13*(n/2))
        if i==2:
            c13 = 0
            for i in range(int(n/2), 2*n-1):
                c13 += vecCharge[i]
                return np.abs(c13*(3/2)*n)
        if i==3:
            c13 = 0
            for i in range(2*n, int(7/2)*n-1):
                c13 += vecCharge[i]
                return np.abs(c13 * (3 / 2)*n)
        if i==4:
            c13 = 0
            for i in range(int(7 / 2) * n, (4*n)-1):
                c13 += vecCharge[i]
                return np.abs(c13*(n/2))

class ErrorNum(BaseException):
    def __init__(self, mensaje):
        print(mensaje)

if __name__ == "__main__":
    t=Tubo(100)
    print(t.calcCoef(3,1))
    #print(t.matrixSearch(3,3))
    #print("Mesh Relaxation method")
    #print("Info:\n\tTamaño de Mesh:",t.gridSize)
    #cij=t.capacitanceCoeff(1,3,0.00001)
    #print("\tC13=\t",cij)
    #l=Tubo(10)
    #cji=l.capacitanceCoeff(3,1,0.00001)
    #print("\tC31=\t",cji)
    #print("\tPromedio(C13+C31/2)= ",(cij+cji)/2)
    #print("\tdiff=\t",np.abs(cij-cji))
