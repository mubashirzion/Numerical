#Root FInding

def bisectionMethod(f, a, b, tol=1e-6, iteration=100):

     if f(a) * f(b) > 0:
        print("Error: f(a) and f(b) must have different signs")
        return None, 0

     c = None        
     prev_c = None     
     iter = 0         

     for i in range(iteration):
        iter += 1

        c = (a + b) / 2 #(Bisection)  b - f(b) * ((b - a) / (f(b) - f(a)) (False, Secant), b - f(b) / df(b)(Newton)
       

        if prev_c == c or abs(f(c)) < tol:
                break

        prev_c = c

        if f(a) * f(c) > 0:
            a = c 
        else:
            b = c 

     return c, iter


#Differ

def backward_method (f,x,h):
    f_backward = f(x-h) # x+h
    f_current = f(x)
    d = (f_current - f_backward)/h #switch
    return d

def f1(x):
    return x**2 -2

x_point = 2
size = 0.01
true_value = 4
d1=backward_method(f1,x_point,size)
print(f" {d1:.10f}")
print(f"{abs(d1  - true_value):.10f}")


#Integ

def Trapoziodal(f,a,b,n):
    h = (b-a)/n
    acc_even =   0
    acc_odd =   0
    for i in range(1,n):
        xi = a+i*h
        if i%2 == 0: # i%3 != 0....not_3 +=f(xi)
            acc_even+=f(xi)
        else:
            acc_odd +=f(xi)
        
    integral =  (h/3) * (f(a)+f(b)+4*acc_odd + 2*acc_even)    #Trap = (h/2) * (f(a) + 2*acc +f(b))  # 3/8 = (3*h / 8) *(f(a)+f(b)+ 3*acc_not_3 + 2 *acc_yes_3)
    return integral

def f1(x):
    return x**2 +1
a,b = 0,2
n=30
exect_value = 4.6667

result  = Trapoziodal(f1,a,b,n)

print(f'{(b-a)/n}')
print(f"{result:.10f}")
print(f"{result - exect_value}")


#Inter

def lagrange(x_list, y_list, x):
    n = len(x_list)
    result = 0

    for i in range(0, n):
        L = 1
        for j in range(0, n):
            if i != j:
                L*= (x-x_list[j])/(x_list[i]-x_list[j])
        
        result += y_list[i]*L

    return result



def newton(x_list, y_list, x):
    n = len(x_list) - 1

    f = []
    for i in range(n+1):
        f.append([0.0]* (n+1))
    
    for i in range(n+1):
        f[i][0] = y_list[i]
    
    for j in range(1, n+1):
        for i in range(n-j+1):
            f[i][j] = (f[i+1][j-1] - f[i][j-1])/(x_list[i+j] - x_list[i])


    p = f[0][0]
    for i in range(1, n+1):
        pd = 1
        for j in range(i):
            pd*= (x-x_list[j])

        p+= f[0][i] * pd

    return p



#Linear
def isStrictlyDiagonallyDominant(A):

    n = len(A)
    dominance = True

    for i in range(n):

        diagonal = abs(A[i][i])

        sum_off_diagonal = 0
        for j in range(n):
            if j != i: 
                sum_off_diagonal += abs(A[i][j])


        if diagonal <= sum_off_diagonal:
            dominance = False

    return dominance


def jacobiMethod(A, b, x0, tol=1e-3, iteration=200):

    if isStrictlyDiagonallyDominant(A):
        print("Will converge")
    else:
        print("Might not converge")

    n = len(A)

    if x0 is None:
        x = [0.0] * n 
    else:
        x = x0.copy() 

    x_new = [0.0] * n
    iter = 0

    for k in range(iteration):

        for i in range(n):

            sum_off_diagonal = 0
            for j in range(n):
                if j != i:
                    sum_off_diagonal += A[i][j] * x[j]

            x_new[i] = (b[i] - sum_off_diagonal) / A[i][i]

        convergence = True
        for i in range(n):
            if abs(x_new[i] - x[i]) > tol:
                convergence = False
                break

        x = x_new.copy()
        iter += 1

        if convergence:
            break

    return x, iter

def gaussSeidelMethod(A, b, x0, tol=1e-3, iterations=200):


    if isStrictlyDiagonallyDominant(A):
        print("Will converge")
    else:
        print("Might not converge")

    n = len(A)

    if x0 is None:
        x = [0.0] * n 
    else:
        x = x0.copy() 

    iter = 0

    for k in range(iterations):

        x_old = x.copy()

        for i in range(n):
            sum_terms = 0

            for j in range(i):
                sum_terms += A[i][j] * x[j] 

            for j in range(i+1, n):
                sum_terms += A[i][j] * x_old[j]

            x[i] = (b[i] - sum_terms) / A[i][i]

        convergence = True
        for i in range(n):
            if abs(x[i] - x_old[i]) > tol:
                convergence = False
                break

        iter += 1

        if convergence:
            break

    return x, iter