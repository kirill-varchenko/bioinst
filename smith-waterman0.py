from pprint import pprint

def smith_waterman0(A, B, mu = -1, delta = -1):
    """
        Smith–Waterman algorithm for local sequience alignment with match score equal +1 and linear gap penalty.

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
    T = [[0] * (n+1) for _ in range(m+1)]

    # Max score and coords
    score_max, i_max, j_max = 0, 0, 0

    # Filling the table
    for i in range(1, n+1):
        for j in range(1, m+1):
            candidates = (T[j-1][i-1] + (1 if A[i-1] == B[j-1] else mu), # go from left-top
                        T[j-1][i] + delta, # go from top
                        T[j][i-1] + delta, # go from left
                        0)
            T[j][i] = max(candidates)
            if T[j][i] > score_max:
                score_max, i_max, j_max = T[j][i], i, j

    # Tracing path back to origin
    i = i_max
    j = j_max
    A_aligned = ''
    B_aligned = ''

    while i > 0 or j > 0:
        if T[j][i] == 0:
            break
        i_min, j_min = i, j
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
    A_global = '-'*(j_min-min(i_min, j_min)) + A[:i_min-1] + A_aligned + A[i_max:] + '-'*(m-j_max)
    B_global = '-'*(i_min-min(i_min, j_min)) + B[:j_min-1] + B_aligned + B[j_max:] + '-'*(n-i_max)
    return A_aligned, B_aligned, A_global, B_global, score_max

A = 'GCATGCU'
B = 'GATTACA'

# A = 'PLEASANTLY'
# B = 'MEANLY'

# A = 'ACTGAG'
# B = 'GCTACT'

A = 'RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGP' 
B = 'CPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFE'

A = 'STIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQECLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVNGVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYA'
B = 'STIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQSTLAQMYPLQEIQNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQECLLLEPGLNEIMANSLDYNERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVNGVDGYDYSRGQLIEDVEHTFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMVDQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINFLLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKSEPWTLALENVVGAKNMNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYAD'

a, b, a_glob, b_glob, t = smith_waterman0(A, B)
print(a_glob)
print(b_glob)
for i in range(len(a_glob)):
    if a_glob[i] != b_glob[i]:
        print(i, a_glob[i], b_glob[i])


