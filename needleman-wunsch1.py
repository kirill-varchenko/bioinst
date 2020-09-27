from pprint import pprint

def needleman_wunsch1(A, B, S, delta = -1):
    """
        Needleman-Wunsch algorithm for sequience alignment with similarity matrix and linear gap penalty.

        Parameters
        ----------
        A, B : str
            Two sequiences
        S : matrix
            Similarity matrix
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
            candidates = (T[j-1][i-1] + S[A[i-1]][B[j-1]], # go from left-top
                        T[j-1][i] + delta, # go from top
                        T[j][i-1] + delta) # go from left
            T[j][i] = max(candidates)

    # Tracing path back to origin
    i = n
    j = m
    A_aligned = ''
    B_aligned = ''

    while i > 0 or j > 0:
        if i > 0 and j > 0 and T[j][i] == T[j-1][i-1] + S[A[i-1]][B[j-1]]:
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

A = 'RKANRK'
B = 'RANNAKA'

S1 = {'A': {'A': 5, 'R': -2, 'N': -1, 'K': -1},
      'R': {'A': -2, 'R': 7, 'N': -1, 'K': 3},
      'N': {'A': -1, 'R': -1, 'N': 7, 'K': 0},
      'K': {'A': -1, 'R': 3, 'N': 0, 'K': 6}}
    
a, b, t = needleman_wunsch1(A, B, S1)
print(t)
print(a)
print(b)

S2 = {'A': {'A': 5, 'R': -2, 'N': -1, 'K': -1},
      'R': {'A': -2, 'R': 7, 'N': -3, 'K': 3},
      'N': {'A': -1, 'R': -3, 'N': 7, 'K': 0},
      'K': {'A': -1, 'R': 3, 'N': 0, 'K': 6}}

a, b, t = needleman_wunsch1(A, B, S2)
print(t)
print(a)
print(b)

pprint(S1)