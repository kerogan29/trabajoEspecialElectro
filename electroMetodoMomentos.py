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
        dist2 = np.power(self.matrixSearch(i1, i2)[1] - self.matrixSearch(int(j1), int(j2))[1], 2)
        dist = np.sqrt(dist1 + dist2)
        if dist==0:
            return 0
        else:
            return np.log(dist)

    def lMatrix(self):
        n=self.gridSize
        matdist= self.matrixDistAux()
        l=np.zeros((2*n,2*n)) # Defino la matriz l como una cuadrada de n/2+n
        ratio=(1/2)*(np.log(1/2))-1/2 # valor usado cuando i=j
        i=0 # Fila
        while i< l.shape[0]:
            j=0 # Columnna
            while j<l.shape[0]:
                if i==j :
                    temp=ratio          #ratio*(1/(s.pi*s.epsilon_0))
                    l[i][j]= temp       # FALTA EPSILON_0
                elif (i<n/2 and j<n/2):
                    l[i][j]= np.log(np.abs(i-j))#1/(s.pi*np.abs(i-j)) # FALTA EPSILON_0
                    l[j][i]=l[i][j]
                elif (i>n/2 and j>n/2):
                    if ((i-n/2)<n and (j-n/2)<n): # Me encuntro en segmento horizontal del conductor 3
                        l[i][j] = np.log(np.abs(i - j))#1 / (s.pi * np.abs(i - j))  # FALTA EPSILON_0
                        l[j][i] = l[i][j]
                    elif ((i-n/2)>n and (j-n/2)>n):
                        l[i][j]= np.log(np.abs(i - j)) #1 / (s.pi * np.abs(i - j))
                        l[j][i]= l[i][j]
                elif (i<n/2 and j>n/2):
                    if j<(3/2)*n:
                        dist1= np.power((matdist[0][i][0][0]-matdist[n-1][j-n/2][0][0]),2)
                        dist2= np.power((matdist[0][i][0][1]-matdist[n-1][j-n/2][0][1]),2)
                        if not (dist1 == 0 and dist2 == 0):
                            l[i][j]= np.log(np.sqrt(dist1+dist2)) #1/(s.pi * np.sqrt(np.power((matdist[0][i][0][0]-matdist[n-1][j-n/2][0][0]),2)+np.power((matdist[0][i][0][1]-matdist[n-1][j-n/2][0][1]),2)))
                            l[j][i]= l[i][j]
                        else:
                            print(matdist[0][i][0][0],matdist[n-1][j-n/2][0][0],matdist[0][i][0][1],matdist[n-1][j-n/2][0][1])
                    else:
                        dist1= np.power((matdist[0][i][0][0] - matdist[j - (3/2)*n][0][0][0]), 2)
                        dist2=  np.power((matdist[0][i][0][1] - matdist[j - (3/2)*n][0][0][1]), 2)
                        if not (dist1== 0 and dist2== 0):
                            l[i][j] = np.log(np.sqrt(dist1+dist2)) #1/(s.pi * np.sqrt(np.power((matdist[0][i][0][0] - matdist[j - (3/2)*n][0][0][0]), 2)+ np.power((matdist[0][i][0][1] - matdist[j - (3/2)*n][0][0][1]), 2)))
                            l[j][i] = l[i][j]
                        else:
                          print(matdist[0][i][0][0],"llol")
                j+=1
            i+=1
        return l

    def llMatrix(self):
        n = self.gridSize
        l = np.zeros((4 * n, 4 * n))  # Defino la matriz l como una cuadrada de n/2+n
        ratio = (1 / 2) * (np.log(1 / 2)) - 1 / 2  # valor usado cuando i=j
        i = 0  # Fila
        while i < l.shape[0]:
            j = i  # Columnna
            while j < l.shape[0]:
                if i == j:
                    temp = ratio  # ratio*(1/(s.pi*s.epsilon_0))
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
                                l[i][j] = self.logDist(i-1,n-1,j-n, n - 1)
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
                            l[i][j] = self.logDist(n-1, i-(5/2)*n, n-1, j-(5/2)*n)
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
class ErrorNum(BaseException):
    def __init__(self, mensaje):
        print(mensaje)

if __name__ == "__main__":
    t=Tubo(4)
    print(t.llMatrix())
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
