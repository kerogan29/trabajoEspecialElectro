''' Trabajo especial de Electrómagnetismo
    Resolución mediante el método de relajación
    Facultad de Ingeniería - Universidad de Montevideo
    2016
'''

# Modules goes here
import numpy as np          # Modulo de calculo matricial
from scipy import constants as s # Modulo para uso de permitividad en el vacio
# Code starts here
class Tubo:                 # Clase que representa a la seccion transversal del tubo
                            # Definir varibles
    def __init__(self, n):  # Constructor; n -> numero de bloques por lado del tubo cuadrado
        self.gridSize=n     # gridsize -> Lado del tubo
        if n % 2 == 0:      # solo puede ser par (tenemos a/2)
            self.seccTubo = np.zeros((n, n)) # Definir seccion tranversal
        else:
            raise ErrorNum("El numero de divisiones debe ser par") # Levanta una excepción si el numero no es par

    def setVoltage(self,nrocond): # Se indica cual conductor estara a V=1v
        n=self.gridSize
        if nrocond==1:      # Cond 1
            for c in range(0,int(n/2)):
                self.seccTubo[0][c]=1
        elif nrocond==2:    # Cond 2
            for c in range(0, int(n/2)-1):
                self.seccTubo[0][int(n/2)+c]= 1
            for c in range(0, n):
                self.seccTubo[c][n-1]= 1
        elif nrocond==3:    # Cond 3
             for c in range(0,int(n/2)-1):
                 self.seccTubo[int(n/2)+c][0]=1
             for c in range(0,n):
                 self.seccTubo[n-1][c]=1
        elif nrocond==4:    # Cond 4
            for c in range(0,int(n/2)):
                self.seccTubo[c][0]=1

    def oneMeshRelaxation(self):
        listaValores=[] # Lista con los valores
        i=0 # Position x
        while i < self.gridSize-1:
            j=0 # Position y
            while j < self.gridSize-1:
                if (i > 0 and j > 0): # Para (i,j) dentro del tubo
                    self.seccTubo[i][j]=(self.seccTubo[i-1][j]+self.seccTubo[i+1][j]+self.seccTubo[i][j-1]+self.seccTubo[i][j+1])/4
                    listaValores.append(self.seccTubo[i][j])
                j+=1
            i+=1
        return listaValores # Devuelvo valores

    def manyMeshRelaxation(self,nro): # Realizo oneMeshRelaxation varias veces
        for i in range(0,nro-1):
             self.oneMeshRelaxation()
        return self.oneMeshRelaxation()

    def errorMeshRelaxation(self, error): # Realiza el procedimiento revisando un error
        count=0
        ctrl=True
        listaOld=[]
        while ctrl:
            listaNew=self.oneMeshRelaxation()
            count+=1
            if listaOld == []:
                listaOld= listaNew
            else:
                ctrl =False
                cond = (np.abs(np.array(listaNew) - np.array(listaOld)) > error) # Controla el error
                listaOld=listaNew
                for i in cond: # Para cada entrada del tubo revisa que sea menor al error
                    if i:
                        ctrl = True # Si hay valores mayores se repite la iteración
        return count

    def chargeDensity(self,nrocond): # Calcula la densidad de carga de un conductor
         ret = np.zeros((self.gridSize, self.gridSize))
         ep0 = s.epsilon_0           # Valor de la permitividad electrica en el vacio
         n=self.gridSize
         if nrocond == 1:            # Valor para el conductor 1
             for c in range(0, int(n / 2)):
                 ret[0][c] = ep0*self.seccTubo[1][c]
         elif nrocond == 2:          # Valor para el conductor 2
             for c in range(0, int(n / 2) - 1):
                 ret[0][int(n / 2) + c] = ep0*self.seccTubo[1][int(n / 2) + c]
             for c in range(0, n):
                 ret[c][n - 1] = ep0*self.seccTubo[c][n - 2]
         elif nrocond == 3:          # Valor para el conductor 3
             for c in range(0, int(n / 2) - 1):
                 ret[int(n / 2) + c][0] = ep0*self.seccTubo[int(n / 2) + c][1]
             for c in range(0, n):
                 ret[n - 1][c] = ep0*self.seccTubo[n - 2][c]
         elif nrocond == 4:          # Valor para el conductor 4
             for c in range(0, int(n / 2)):
                 ret[c][0] = ep0*self.seccTubo[c][1]
         return (ret,nrocond) # Devuelve tupla con la matriz de la densidad de carga y el nro de cond

    def capacitanceCoeff(self, i, j, error): # Coeficiente de capacitancia Cij
        nrocond= j  # Nro del conductor
        sum1= 0      # Suma de las densidades de carga
        sum2= 0      # SUma de las densidades de carga
        coeff= 0    # Coefciente
        n= self.gridSize    # Tamaño del tubo
        self.setVoltage(i)
        self.errorMeshRelaxation(error)
        ret= self.chargeDensity(j)[0]
        # Todos los conductores tiene una carga igual a la suma de las densidades lineales
        # multiplicadas por la longitud del conductor y dividida entre un potencial de 1v
        if nrocond == 1:
            for c in range(0, int(n / 2)):
                sum1+= ret[0][c]
            coeff= (sum1*(n/2))
        elif nrocond == 2:
            for c in range(0, int(n / 2) - 1):
                sum1 += ret[0][int(n / 2) + c]
            coeff = (sum1 * (n / 2))
            for c in range(0, n):
                sum2 += ret[c][n - 1]
            coeff += (sum2 * n)
        elif nrocond == 3:
            for c in range(0, int(n / 2) - 1):
                sum1 += ret[int(n / 2) + c][0]
            coeff = (sum1 * (n / 2))
            for c in range(0, n):
                sum2 += ret[n - 1][c]
            coeff += (sum2 * n)
        elif nrocond == 4:
            for c in range(0, int(n / 2)):
                sum1 += ret[c][0]
            coeff = (sum1 * (n / 2))
        return coeff

class ErrorNum(BaseException):
    def __init__(self, mensaje):
        print(mensaje)

if __name__ == "__main__":
    t=Tubo(10)
    print("Mesh Relaxation method")
    print("Info:\n\tTamaño de Mesh:",t.gridSize)
    cij=t.capacitanceCoeff(1,3,0.0000000000000000000000000000000000000000000001)
    print("\tC13=\t",cij)
    l=Tubo(10)
    cji=l.capacitanceCoeff(3,1,0.0000000000000000000000000000000000000000000001)
    print("\tC31=\t",cji)
    print("\tPromedio(C13+C31/2)= ",(cij+cji)/2)
    print("\tdiff=\t",np.abs(cij-cji))
    #print(t.seccTubo)
    #t.setVoltage(4)
    #t.manyMeshRelaxation(2)
    #print(t.chargeDensity(3)[0])
    #t.setVoltage(1)
    #print(t.seccTubo)






