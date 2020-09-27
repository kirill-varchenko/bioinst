from pprint import pprint

def needleman_wunsch0(A, B, mu = -1, delta = -1):
    """
        Needleman-Wunsch algorithm for sequience alignment with match score equal +1.

        Parameters
        ----------
        A, B : str
            Two sequiences
        mu : num, optional
            Mismatch penalty (default is -1)
        delta : num, optional
            Gap penalty (default is -1)
        """

    n = len(A)
    m = len(B)

    # Initial score table
    T = [[delta*i for i in range(n+1)] if j == 0 else [delta*j] + [0]*n for j in range(m+1)]

    # Filling the table
    for i in range(1, n+1):
        for j in range(1, m+1):
            candidates = (T[j-1][i-1] + (1 if A[i-1] == B[j-1] else mu), # go from left-top
                        T[j-1][i] + delta, # go from top
                        T[j][i-1] + delta) # go from left
            T[j][i] = max(candidates)


    # Tracing path back to origin
    i = n
    j = m
    A_aligned = ''
    B_aligned = ''

    while i > 0 or j > 0:
        if i > 0 and j > 0 and T[j][i] == T[j-1][i-1] + (1 if A[i-1] == B[j-1] else mu):
            # Came from left-top 
            A_aligned = A[i-1] + A_aligned
            B_aligned = B[j-1] + B_aligned
            i -= 1
            j -= 1 
        elif j > 0 and T[j][i] == T[j-1][i] + delta:
            # Came from top 
            A_aligned = '_' + A_aligned
            B_aligned = B[j-1] + B_aligned
            j -= 1
        else:
            # Came from left
            A_aligned = A[i-1] + A_aligned
            B_aligned = '_' + B_aligned
            i -= 1

    # Return two aligned sequiences and total score
    return A_aligned, B_aligned, T[m][n]

A = 'GCATGCU'
B = 'GATTACA'

A = 'PLEASANTLY'
B = 'MEANLY'

#A, B = input().split(' ')

a, b, t = needleman_wunsch0(A, B)
print(a, b)


