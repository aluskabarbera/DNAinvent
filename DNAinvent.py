import colorama

colorama.init()

print(
    """
            DNAinvent       v1.1
        [0]
            1) DNA sequence and his double helix.
            2) Transcription DNA to RNA.
            3) Translate RNA to amino acids.
            4) CRISPR/cas9
            5) Convert DNA to binary.
        [*]

    """
)

option = input(str("Please select a option: "))

if option == "1":
#  DNA sequence and his double helix.
    def doublehelix(dna_string: str) -> str:
        complementary = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
        output = ''
        for nucleotid in dna_string:
            output += complementary[nucleotid]
        return output

    print(doublehelix(input("Input your DNA chain (A, C, G, T): ").upper()))

elif option == "2":
#  Transcription DNA to RNA.
    def convertdnatorna(dna_string: str) -> str:
        RNA = {
        'A': 'A',
        'T': 'U',
        'C': 'C',
        'G': 'G'
    }
        output = ''
        for nucleotid in dna_string:
            output += RNA[nucleotid]
        return output

    print(convertdnatorna(input("Input your DNA chain (A, C, G, T): ").upper()))

elif option == "3":
#  Translate RNA to aminoacids.
    import collections

    aminoacids = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
    }
    key = {
    'A':'Alanine (Ala)',    'R':'Arginine (Arg)',       'N':'Asparagine (Asn)',
    'D':'Aspartate (Asp)',  'C':'Cysteine (Cys)',       'Q':'Glutamine (Gin)',
    'E':'Glutamate (Glu)',  'G':'Glycine (Gly)',        'H':'Histidine (His)', 
    'I':'Isoleucine (Ile)', 'L':'Leucine (Leu)',        'K':'Lysine (Lys)', 
    'M':'Methionine (Met)', 'F':'Phenylalanine (Phe)',  'P':'Proline (Pro)',
    'S':'Serine (Ser)',     'T':'Threonine (Thr)',      'W':'Tryptophan (Trp)', 
    'Y':'Tyrosine (Tyr)',   'V':'Valine (Val)'
    }

    def get_input():
        return(input(
            "Input your RNA chain (A, C, G, U): "
        ).upper())
        
#· For 3 amino acid in DNA sequence, we obtain one codon.
    def cutted_surplus_translated(DNA):
        tmp = DNA[:len(DNA) - len(DNA)%3]
        protein = [aminoacids.get(tmp[i:i+3]) for i in range(0, len(tmp), 3)]
        return(protein)

#· This codon coincide with any protein?. 
    def print_frequencies(frequencies):
        for k, v in frequencies.items():
            print(f"The protein {k} is repeated {v} times")

#· Scientific name of codons
    def get_proteins():
        DNA = get_input()
    
        protein = cutted_surplus_translated(DNA)
        print(protein)
    
        print_frequencies(collections.Counter(protein))
    
        return([key.get(i) for i in protein])

    print(get_proteins())

elif option == "4":
#  CRISPR/cas9
    DNA = input("Input your DNA chain (A, C, G, T): ")

#· Write the sequence of nitrogenous bases you want to change and write the new one. 
    def cut(DNA_chain):
        cutter = input("Insert the sequence that you want cut: ")
        if cutter not in DNA:
            print("This sequence is not in the chain.")
        change = input("Insert the sequence that you want paste: ")
        return DNA_chain.replace(cutter, f'\033[96m{change}\033[0m')

    new_chain = cut(DNA)
    print(new_chain)

elif option == "5":
#· Convert DNA to binary.

    def convertdnatobinary(dna_string: str) -> str:
        d = {'A':'0001','T':'0010','G':'0011','C':'0100',}
        output = ''
        for nucleotid in dna_string:
            output += d[nucleotid]
            if nucleotid not in d:
                return(nucleotid + " is not a valid character. Run this option again.")
        return output

    print(convertdnatobinary(input("Input your DNA chain (A, C, G, T): ").upper()))

elif option == "*":
    print(
    """
            DNAtools       v1.0
        [0]
            1) Find DNA sequence in the original sequence.
            2) Compare whether the first or second amino acid sequence is: greater, lesser, equal, greater than or equal, lesser than or equal.
            3) Know the percentage that two DNA sequences have in common.
            4) Generate permutations from the DNA sequence.

    """
)

    option2 = input(str("Please select a option: "))

    if option2 == "1":
    #· Find DNA sequence in the original sequence.
        str1 = input("Input your DNA chain (A, C, G, T): ")
        str2 = input("Insert the sequence you want to find: ")
        str1 = str1.replace(str2,f'\033[96m{str2}\033[0m')
        print(str1)

    elif option2 == "2":
    #· Compare whether the first or second amino acid sequence is: greater, lesser, equal, greater than or equal, lesser than or equal.
        DNA1 = input("Input the first DNA chain (A, C, G, T): ")
        DNA2 = input("Input the second DNA chain (A, C, G, T): ")

        print('The first DNA sequence is greater than the second DNA sequence:', DNA1.upper() > DNA2.upper())
        print('The first DNA sequence is lesser than the second DNA sequence:', DNA1.upper() < DNA2.upper())
        print('The first DNA sequence is equal to the second DNA sequence:', DNA1.upper() == DNA2.upper())
        print('The first DNA sequence is greater than or equal to the second DNA sequence:', DNA1.upper() >= DNA2.upper())
        print('The first DNA sequence is lesser than or equal to the second DNA sequence:', DNA1.upper() <= DNA2.upper())
    
    elif option2 == "3":
    #· This function returns the match percentage of two DNA sequences.

        first = str(input("Input your dna sequence (A, G, C, T): "))
        second = str(input("Input other dna sequence to campare with another sequence: "))

        def similar(first, second):
            accuracy = 0
            while len(first) > len(second):
                second = second + " "
            while len(second) > len(first):
                first = first + " "
            for i, l in enumerate(first):
                if second[i] == l:
                    accuracy += 1
            return accuracy/((len(first)+len(second))/2)*100

        print(similar(first, second), "%")
    
    elif option2 == "4":
        # Permutations.
        from itertools import permutations
        DNA = input("Input your dna sequence (A, G, C, T): ") 

        letras = [letra for letra in DNA if letra in 'ACTG' or 'actg']
        print(letras)
        print()

        permutaciones = [p for p in permutations(letras, 2)]

        contador = 0

        for p in permutaciones:
            print(p)
            contador += 1
        print()

        print(DNA)
        print(f'Number of permutations generated: {contador}.')
