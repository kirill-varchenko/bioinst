from pprint import pprint

def needleman_wunsch2(A, B, kappa = 1, mu = -1, delta = -1, epsilon = -1):
    """
        Needleman-Wunsch algorithm for sequience alignment with affine gap penalty (Gotoh method).

        Parameters
        ----------
        A, B : str
            Two sequiences
        kappa : num, optional 
            Match reward (default is +1)
        mu : num, optional
            Mismatch penalty (default is -1)
        delta : num, optional
            Gap open penalty (default is -1)
        epsilon : num, optional
            Gap extend penalty (default is -1)
        """

    n = len(A)
    m = len(B)

    # Best score table
    D = [[delta + (i-1)*epsilon for i in range(n+1)] if j == 0 else [delta + (j-1)*epsilon] + [0]*n for j in range(m+1)]
    D[0][0] = 0
    # Score if x_i aligns to a gap after y_j
    P = [[-1000] + [0] * (n) for _ in range(m+1)]
    P[0][0] = 0
    # Score if y_j aligns to a gap after x_i
    Q = [([-1000] * (n+1) if j == 0 else [0] * (n+1)) for j in range(m+1)]
    Q[0][0] = 0
 
    # Filling the table
    for i in range(1, n+1):
        for j in range(1, m+1):
            P[j][i] = max(D[j][i-1]+delta, P[j][i-1]+epsilon)
            Q[j][i] = max(D[j-1][i]+delta, Q[j-1][i]+epsilon)
            D[j][i] = max(D[j-1][i-1] + (kappa if A[i-1] == B[j-1] else mu), 
                        P[j][i], 
                        Q[j][i])

    # Tracing path back to origin
    i = n
    j = m
    A_aligned = ''
    B_aligned = ''
    table = 'D'

    while i > 0 or j > 0:
        if i == 0:
            table = 'Q'
        elif j == 0:
            table = 'P'

        if table == 'D':
            if D[j][i] == Q[j][i]:
                # Switch table
                table = 'Q'
            elif D[j][i] == P[j][i]:
                # Switch table
                table = 'P'
            elif D[j][i] == D[j-1][i-1] + (kappa if A[i-1] == B[j-1] else mu):
                # Came from left-top 
                A_aligned = A[i-1] + A_aligned
                B_aligned = B[j-1] + B_aligned
                i -= 1
                j -= 1
            
        elif table == 'P':
            if P[j][i] == D[j][i-1] + delta:
                table = 'D'
            A_aligned = A[i-1] + A_aligned
            B_aligned = '_' + B_aligned
            i -= 1
        elif table == 'Q':
            if Q[j][i] == D[j-1][i] + delta:
                table = 'D'
            A_aligned = '_' + A_aligned
            B_aligned = B[j-1] + B_aligned
            j -= 1

    # Return two aligned sequiences and total score
    return A_aligned, B_aligned, D[m][n]

# A = 'GCATGCU'
# B = 'GATTACA'

A = 'PLEASANTLY'
B = 'MEANLY'

#A = 'ATGTAGTGTATAGTACATGCA'
#B = 'ATGTAGTACATGCA'

#A = 'GAAAAAAT'
#B = 'GAAT'

A = 'ACCT'
B = 'CC'


a, b, t = needleman_wunsch2(A, B, mu=-1, delta=-3, epsilon=-1)
print(a)
print(b)


