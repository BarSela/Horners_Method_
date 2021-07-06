import sympy as sp
import numpy as np
from sympy.utilities.lambdify import lambdify
x = sp.symbols('x')

def polynomial(p):
    """
    The function get a list of Coefficients and create Polynom
    :param p: list of Coefficients
    :return:polynom
    """
    n=len(p)-1
    poly=0
    for i in range(0,len(p)):
        if(n>=2):
            poly+=p[i]*x**n
        if(n==1):
            poly += p[i]*x
        if (n == 0):
            poly += p[i]
        n-=1
    return poly

def find_root(p,q, guess ):
    """
    The function find roots by approximation to the root.
    x_n it's the first guess
    :param p: polynom
    :param n: degree of the polynom
    :return: root of p
    """
    p_calc = lambdify(x, p)
    q_calc = lambdify(x, q)
    x_n = guess - (p_calc(guess) / q_calc(guess))

    return x_n


def div_poly(p,guess):
    """
    synthetic division of p in x-root to simplify the polynom
    :param p: the origin polynom
    :param root:the root of p that was found
    :return:new polynom (degree n-1)
    """
    x = sp.symbols('x')
    poly = sp.poly(p).all_coeffs()
    result =[]
    result.append(poly[0])
    for i in range(1,len(poly)-1):
        x=poly[i]+result[i-1]*guess
        result.append(x)

    return polynomial(result)



def Horners_Method(p,n,guess):
    """
    find roots of polynom
    the function print to the screen the roots of the polynon (If there are any at all)
    :param p: the polynom
    :param n: the degree of the polynom
    """
    result = []  # contain all the roots that was found
    e = 0.000001  # The required level of accuracy
    limit=[-1000,1000]  # Setting boundaries for testing in case of not finding a root
    p_calc = lambdify(x, p)

    for i in range(0,n):
        root = guess
        flag=True
        j=1
        print('Try Finding a root by approximation:')
        while abs(p_calc(root)) > e:
            Q_x = div_poly(p, root)
            root=find_root(p,Q_x,root)
            print('approximation ',j,': ', root)
            if root<limit[0] or root>limit[1] or j>50:
                flag=False
                print('No more roots ')
                break
            j += 1

        p=Q_x
        p_calc = lambdify(x, p)

        if flag:
            result.append(round(root, 5))  # round the result according to the selected epsilon
            print("found root: ", root)
        else:
            break
    if len(result)==0:
        print("No roots")
    else:
        print(result)
        return result



f = x**5-8*x**4-72*x**3+382*x**2+727*x-2310
f1= 2*x**4-3*x**2+3*x-4
f2=x**2-6*x+8
f3=5*x**4-2*x**2+7*x-4
f4=x**2+4*x+6
f5=5*x**5+4*x**3-2*x
Horners_Method(f2,2,0)
