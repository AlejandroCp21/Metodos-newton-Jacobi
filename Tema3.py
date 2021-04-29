#!/usr/bin/env python
# coding: utf-8

# ## Método de Jacobi
# 
# ${x_i^k} = T*x_{k-1} + c$
# 

# In[9]:


import numpy as np

from math import *


# In[30]:


#Definimos la matriz T
T = np.array([[0,1/10,-1/5,0],[1/11,0,1/11,-3/11],[-1/5,1/10,0,1/10],[0,-3/8,1/8,0]])
c = np.array([3/5,25/11,-11/10,15/8])
x = np.array([0,0,0,0]) #vector de valores iniciales
erroraceptado=0.01


# In[31]:


print("Matriz de :", T.shape)
print(c.shape)
error=1
while error>0.01:
    resultado=np.dot(T,x) + c
    print(resultado)
    numerador = np.amax( abs((resultado - x)))
    denominador = np.amax( abs(resultado))
    error = numerador/denominador
    x=resultado
    


# In[ ]:





# In[25]:


numerador = np.amax( abs((resultado - x)))
denominador = np.amax( abs(resultado))
error = numerador/denominador
print(error)


# ## Método de Newton
# 

# In[46]:


def jacobiano(x):
    J = np.array([[3,x[2]*sin(x[1]*x[2]),x[1]*sin(x[1]*x[2])],[2*x[0], -162*(x[1]+0.1), cos(x[2])],[-x[1]*exp(-x[0]*x[1]), -x[0]*exp(-x[0]*x[1]), 20]])
    JInversa = np.linalg.inv(J)
    return JInversa

def Fx(x):
    xk=np.array([3*x[0]-cos(x[1]*x[2])-1/2, x[0]**2-81*(x[1]+0.1)**2+sin(x[2])+1.06, exp(-x[0]*x[1])+ 20*x[2] + (10-math.pi-3)/3])
    
    #x1 = 3*x[0]-cos(x[1]*x[2])-1/2
    
    #x2= x[0]**2-81*(x[1]+0.1)**2+sin(x[2])+1.06
    
    #x3=exp(-x[0]*x[1])+ 20*x[2] + (10-math.pi-3)/3
    #xk = np.array([x1,x2,x3])
    return xk


# In[59]:


x = np.array([0.1,0.1,-0.1]) #vector de valores iniciales
error = 1
c=0


# In[60]:


while error>0.01:
    c+=1
    r= jacobiano(x)
    feval = Fx(x)
    resultado =x- np.dot(r,feval)
    
    numerador = np.amax( abs((resultado - x)))
    denominador = np.amax( abs(resultado))
    error = numerador/denominador
    x=resultado
    
    print("Iteración:",c," Resultado ",resultado)


# In[ ]:




