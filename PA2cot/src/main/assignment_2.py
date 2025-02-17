import numpy as np  

def neville_method():
    x = [3.6, 3.8, 3.9]
    val = [1.675, 1.436, 1.318]
    w = 3.7  

    n = len(x)
    neville = np.zeros((n, n))
    neville[:, 0] = val

    for i in range(n):
        neville[i][0] = val[i]

    for i in range(1, n):
        for j in range(1, i + 1):
            term1 = (w - x[i - j]) * neville[i][j - 1]
            term2 = (w - x[i]) * neville[i - 1][j - 1]
            neville[i][j] = (term1 - term2) / (x[i] - x[i - j])

    for i in range(n):
        row_str = ""
        for j in range(i + 1):
            row_str += f"{neville[i][j]:0.16f} "
        print(row_str)


def newton_forward():
    xi = np.array([7.2, 7.4, 7.5, 7.6])
    fxi = np.array([23.5492, 25.3913, 26.8224, 27.4589])

    lim = len(xi)
    diffs = np.zeros((lim, lim))

    diffs[:, 0] = fxi

    for i in range(1, lim):
        for j in range(1, i+1):
            diffs[i, j] = (diffs[i, j-1] - diffs[i-1, j-1]) / (xi[i] - xi[i - j])

    print("Divided difference table:\n")
    for i in range(lim):
        row_str = ""
        for j in range(i + 1):
            if diffs[i, j] >= 0:
                row_str += " "
            row_str += f"{diffs[i, j]:0.16f} "
        print(row_str)

    x_target = 7.3
    approx = diffs[0, 0]  
    product = 1.0
    for j in range(1, lim):
        product *= (x_target - xi[j-1])
        approx += product * diffs[j, j]

    print(f"\nApproximation f({x_target}) = {approx:0.16f}")


def newton_formula(x, f, x_target):
    n = len(x)
    table = np.zeros((n, n))
    table[:, 0] = f

    for j in range(1, n):
        for i in range(n - j):
            table[i, j] = (table[i+1, j-1] - table[i, j-1]) / (x[i+j] - x[i])

    result = table[0, 0]
    product = 1.0
    for i in range(1, n):
        product *= (x_target - x[i-1])
        result += product * table[0, i]

    return result


def Hermite_poly_matrix():
    np.set_printoptions(
    formatter={'float': '{: .8e}'.format},  # 8 decimal places, exponent form
    floatmode='maxprec_equal',
    linewidth=120
    )

    x = np.array([3.6, 3.8, 3.9])
    fx = np.array([1.675, 1.436, 1.318])
    dfx = np.array([-1.195, -1.188, -1.182])

    n = len(x)
    zval = np.zeros(2 * n)
    for i in range(n):
        zval[2 * i] = x[i]
        zval[2 * i + 1] = x[i]

    #print("zval array:")
    #print(zval)

    fz = np.zeros(2 * n)
    for i in range(n):
        fz[2 * i] = fx[i]
        fz[2 * i + 1] = fx[i]

    #print("fz array:")
    #print(fz)

    fdd = np.zeros(2*n)
    sdd= np.zeros(2*n)
    tdd = np.zeros(2*n)

    fdd[0] = 0.0
    fdd[1] = dfx[0]

    #fill in first dd second occurence of node is deriv
    fdd[1] = dfx[0]   
    fdd[3] = dfx[1]  
    fdd[5] = dfx[2]  

    # Compute standard divided differences between distinct nodes
    fdd[0] = 0.0  
    fdd[2] = (fz[2] - fz[1]) / (zval[2] - zval[1])  
    fdd[4] = (fz[4] - fz[3]) / (zval[4] - zval[3])  

    #print("fdd array:")
    #print(fdd)

    tol = 1e-14

    sdd = np.zeros(2*n)
    sdd[0] = 0.0
    sdd[1] = 0.0


    sdd[2] = (fdd[2] - fdd[1]) / (zval[2] - zval[0])
    sdd[3] = (fdd[3] - fdd[2]) / (zval[3] - zval[1])
    sdd[4] = (fdd[4] - fdd[3]) / (zval[4] - zval[2])
    sdd[5] = (fdd[5] - fdd[4]) / (zval[5] - zval[3])

    #print("sdd array:")
    #print(sdd)

    tdd = np.zeros(2*n)
    tdd[0] = 0.0
    tdd[1] = 0.0
    tdd[2] = 0.0
    tdd[3] = (sdd[3] - sdd[2]) / (zval[3] - zval[0])
    tdd[4] = (sdd[4] - sdd[3]) / (zval[4] - zval[1])
    tdd[5] = (sdd[5] - sdd[4]) / (zval[5] - zval[2])
    #print("tdd array:")
    #print(tdd)

    H_matrix = np.zeros((6, 5))

    for i in range(6):
        H_matrix[i][0] = zval[i]
        H_matrix[i][1] = fz[i]
        H_matrix[i][2] = fdd[i]
        H_matrix[i][3] = sdd[i]
        H_matrix[i][4] = tdd[i]

    print(H_matrix)
def cubic_spline():
    np.set_printoptions(formatter={'float': '{: .8f}'.format}, suppress=True)
    x = np.array([2, 5, 8, 10])
    fx = np.array([3, 5, 7, 9])
    n = len(x)
    h = [x[i+1] - x[i] for i in range(n-1)]
    
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    A[0, 0] = 1 
    A[-1, -1] = 1 
    
    for i in range(1, n-1):
        h_prev = h[i-1]
        h_curr = h[i]
        A[i, i-1] = h_prev
        A[i, i] = 2 * (h_prev + h_curr)
        A[i, i+1] = h_curr
        term1 = (fx[i+1] - fx[i]) / h_curr
        term2 = (fx[i] - fx[i-1]) / h_prev
        b[i] = 3 * (term1 - term2)
    
    c = np.linalg.solve(A, b)
    
    print("Matrix A:")
    print(A)
    print("\nVector b:")
    print(b)
    print("\nVector x (c):")
    print(c)


