from numpy import linalg as LA
from time import process_time
import numpy as np
import secrets
from time import process_time



#******************************************************************
#************BASIC FUNCTION MAX-PLUS ALGEBRA***********************
#******************************************************************

#generate random matrix

def generate_random_matrix(m,n, min, max):
    H = [[secrets.randbelow((max - min) + 1) + min for i in range(m)] for j in range(n)]
    return H

#generate random exponent
def generate_exponent(order):
  secretsGenerator = secrets.SystemRandom()
  m = secretsGenerator.randint(3,2**order)
  return m

#Tropical Identity (Addition)/have checked(correct)
def maxpluszeros(n):
  return np.zeros((n,n))+-np.inf

#Tropical Identity (Multipication)/have checked(correct)
def maxeins(n):
  #id=np.zeros((n,n))
  id=generate_random_matrix(n,n,0,0)
  for i in range(0,n):
    for j in range(0,n):
      if i==j:
        id[i][j]=id[i][j]
      else:
        id[i][j]=-np.inf
  return id
#Tropical addition/have checked it(correct)
#def oplus(A,B):
  #T=np.maximum(A,B)
 # return T  
# Tropical Matrix Addition
def oplus(x, y):
    return [[a if a > b else b for a, b in zip(x_rows, y_rows)] for x_rows, y_rows in zip(x, y)]

def sub(x, y):
    return [[1 if (abs(a-b) <=0.0000001) else 0 for a, b in zip(x_rows, y_rows)] for x_rows, y_rows in zip(x, y)]

def minus(p,q):
  return[[a-b for a,b  in zip(p_rows, q_rows)] for p_rows, q_rows in zip(p, q)]

def plus(p,q):
  return[[a+b for a,b  in zip(p_rows, q_rows)] for p_rows, q_rows in zip(p, q)]

def multi(p,q):
  return[[a*b for a,b  in zip(p_rows, q_rows)] for p_rows, q_rows in zip(p, q)]
    
# Tropical Matrix Multiplication
def otimes(x, y):
    return [[max(a + b for a, b in zip(row, col)) for col in zip(*y)] for row in x]

# Tropical Adjoint Multiplication
def adj_multiply(x, y):
    return oplus(oplus(x, y), otimes(x, y))


# Tropical Adjoint Multiplication
def adj_multiply(x, y):
    return oplus(oplus(x, y), otimes(x, y))


# First semigroup operation from: Tropical cryptography II: Extensions by homomorphisms
def exoticproduct(x, g, y, h):
    return oplus(adj_multiply(x, h), y), adj_multiply(g, h)


# Computes intermediaries with square-and-multiply method utilising associative property of semigroup operation
def adjointpowerfast(M, H, a):
    S = [(M, H)]
    for i in range(len(bin(a)[2:]) - 1):
        S.append(exoticproduct(*S[i], *S[i]))
    temp = True
    result = None
    for i in range(len(bin(a)[2:])):
        if bin(a)[:1:-1][i] == "1" and temp:
            result = S[i]
            temp = False
        elif bin(a)[:1:-1][i] == "1":
            result= exoticproduct(*result, *S[i])
    return result


#faster max plus multiplication
def fastpower(A,t):
  exp = bin(t)
  C = A
  for i in range(3, len(exp)):
    C = otimes(C,C)
    if (exp[i:i+1]=='1'):
      C = otimes(C,A)
  return C


#column 
def column(B,z):
  n=len(B)
  #C=np.zeros((n,n))
  C=generate_random_matrix(n,n,0,0)
  for i in range(n):
    if z[i]==1:
      for k in range(n):
          C[k][i]=B[k][i]
    else:
      for k in range(n):
          C[k][i]=-np.inf
  return C

#row 
def row(B,z):
  n=len(B)
  R=np.zeros((n,n))
  for i in range(n):
    if z[i]==1:
        R[i,:]=B[i,:]
    else:
        R[i,:]=-np.inf
  return R

#maxpower (slow version)
def maxpower(A,t):
  n=len(A)
  C=maxeins(n)
  for i in range(t):
    C=otimes(A,C)
  return C

# detect cycles
def cycles(A,eigvector): 
  n=len(A)
  i=0

  #J=np.zeros((1,n+1))
  J=(n+1)*[0]
  #print('J=')
  for j in range(0,n+1):
    C=A[i,:]+np.transpose(eigvector)
    
    index =np.where(C == max(C))
    #print('index=',index)
    ind=index[0][0]
    i=ind
    #index=mlab.find(C==max(C),1)
    #index =find(C==max(C),1)
    #i=index+np.ones((1,1))
    #print(i)
    J[j]=i+1

  return J
def cyclepolicy(pi):
  n=len(pi)
  J=(n+1)*[0]
  index=0 
  newindex=0 
  #newindex=pi[index]
  #print('newindex',newindex)    
  for j in range(0,n+1):
    index=newindex
    #print('index',index)
    newindex=pi[index]
    J[j]=newindex+1  
  return J
    

def critical(J):
  print('J=',J)
  n=len(J)
  #print(n)
  temp=0
  for j in range(1,n):
    for i in range(j):
      if J[i]==J[j]:
        B=J[i:j]
        print(B)
        temp=1
        break
    if temp==1:
      break
  return B

def maxfloyd(A):
  n=len(A)
  eps=0.0000000001
 
  g=A
  for k in range(n):
     for i in range(n):
       for j in range(n):
         if (i!=k) and (j!=k):
           if g[i][k]+g[k][j]>g[i][j]:
            g[i][j]=g[i][k]+g[k][j]
         if i==j and g[i][i]>eps:
           #print('g[ii]=',g[i][i])
           raise SyntaxError('Terminated after detecting a positive cycle')
  return g


def generate_integer(below,upper):
  secretsGenerator = secrets.SystemRandom()
  m = secretsGenerator.randint(below,upper)
  return m

def diagonalmax(n):
  A=maxeins(n)
  n=len(A)
  for i in range(n):
    for j in range(n):
      A[i][i]=(-1)*generate_integer(1,100)+generate_integer(1,100)
  return A
