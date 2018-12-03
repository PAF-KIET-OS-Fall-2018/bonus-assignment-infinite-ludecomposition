import pprint
import time
from threading import*
PA=None
P=None
L=None
U=None
def mult_matrix(M, N):
    """Multiply square matrices of same dimension M and N"""

    # Converts N into a list of tuples of columns
    tuple_N = zip(*N)

    # Nested list comprehension to calculate matrix multiplication
    NX= [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in N] for row_m in M]
    global PA
    PA=NX

def pivot_matrix(M):
    """Returns the pivoting matrix for M, used in Doolittle's method."""
    m = len(M)

    # Create an identity matrix, with floating point values
    id_mat = [[float(i ==j) for i in range(m)] for j in range(m)]

    # Rearrange the identity matrix such that the largest element of
    # each column of M is placed on the diagonal of of M
    for j in range(m):
        row = max(range(j, m), key=lambda i: abs(M[i][j]))
        if j != row:
            # Swap the rows
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat
def lu_value(start, n,lenth) :
    global L
    global U
    global P
    global PA
    
    # Perform the LU Decomposition
    for j in range(start,n):
        # All diagonal entries of L are set to unity
        L[j][j] = 1.0

        # LaTeX: u_{ij} = a_{ij} - \sum_{k=1}^{i-1} u_{kj} l_{ik}
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = PA[i][j] - s1

        # LaTeX: l_{ij} = \frac{1}{u_{jj}} (a_{ij} - \sum_{k=1}^{j-1} u_{kj} l_{ik} )
        for i in range(j,lenth):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]    


def lu_decomposition(A):
    """Performs an LU Decomposition of A (which must be square)
    into PA = LU. The function returns P, L and U."""
    global L
    global U
    global P
    global PA    
    n = len(A)

    # Create zero matrices for L and U
    L = [[0.0] * n for i in range(n)]
    U = [[0.0] * n for i in range(n)]

    # Create the pivot matrix P and the multipled matrix PA
    P = pivot_matrix(A)
    #PA = mult_matrix(P, A)
    obj=Thread(target=mult_matrix,args=(P,A))
    
    obj.start()
    obj.join;
    threads = []
    for i in range(1,n+1):
        thread = Thread(target=lu_value,args=(i-1,i,n))
        threads.append(thread)
        thread.start()
    for thread in threads:  
        thread.join()       
    
    

    return (P, L, U)



start_time = time.time()

A = [[7, 3, -1, 2], [3, 8, 1, -4], [-1, 1, 4, -1], [2, -4, -1, 6]]
P1, L, U = lu_decomposition(A)
print("--- %s seconds ---" % (time.time() - start_time))
print ("A:")
pprint.pprint(A)

print ("P:")
pprint.pprint(P1)
print ("PA:")
pprint.pprint(PA)

print ("L:")
pprint.pprint(L)

print ("U:")
pprint.pprint(U)