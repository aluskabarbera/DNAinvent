import colorama
import random

colorama.init()

while True:
    print(
    """
            DNAinvent
        [0]
            1) DNA sequence and his double helix.
            2) Transcription DNA to RNA.
            3) Translate RNA to amino acids.
            4) CRISPR/cas9
            5) Convert DNA to binary.
        [*]

    """)
    option = input(str("Please select a option: "))

    if option == "1":
    #  DNA sequence and his double helix.
        def doublehelix(dna_string: str) -> str:
            characters_str = str(len(dna_string))
            print("This sequence have " + characters_str + " characters")
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

        print(doublehelix(input("Input the DNA chain (A, C, G, T): ").upper()))

    elif option == "2":
    #  Transcription DNA to RNA.
        def convertdnatorna(dna_string: str) -> str:
            characters_str = str(len(dna_string))
            print("This sequence have " + characters_str + " characters")
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

        print(convertdnatorna(input("Input the DNA chain (A, C, G, T): ").upper()))

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
                "Input the RNA chain (A, C, G, U): "
            ).upper())
        characters_str = str(len(get_input()))
        print("This sequence have " + characters_str + " characters")
            
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
            characters_str = str(len(DNA))
            print("This sequence have " + characters_str + " characters")
            protein = cutted_surplus_translated(DNA)
            print(protein)
            print_frequencies(collections.Counter(protein))
            return([key.get(i) for i in protein])
        print(get_proteins())

    elif option == "4":
    #  CRISPR/cas9
        DNA = input("Input the DNA chain (A, C, G, T): ")
        characters_str = str(len(DNA))
        print("This sequence have " + characters_str + " characters")


    #· Write the sequence of nitrogenous bases you want to change and write the new one. 
        def cut(DNA_chain):
            cutter = input("Insert the sequence that you want cut: ")
            characters_str_cutter = str(len(cutter))
            print("This sequence have " + characters_str_cutter + " characters")
            if cutter not in DNA:
                print("This sequence is not in the chain.")
                exit()
            change = input("Insert the sequence that you want paste: ")
            characters_str_change = str(len(change))
            print("This sequence have " + characters_str_change + " characters")
            return DNA_chain.replace(cutter, f'\033[96m{change}\033[0m')
        new_chain = cut(DNA)
        print(new_chain)

    elif option == "5":
    #· Convert DNA to binary.
        def convertdnatobinary(dna_string: str) -> str:
            characters_str = str(len(dna_string))
            print("This sequence have " + characters_str + " characters")
            d = {'A':'0001','T':'0010','G':'0011','C':'0100',}
            output = ''
            for nucleotid in dna_string:
                output += d[nucleotid]
                if nucleotid not in d:
                    return(nucleotid + " is not a valid character. Run this option again.")
            return output
        print(convertdnatobinary(input("Input the DNA chain (A, C, G, T): ").upper()))

    elif option == "*":
        print(
        """
                DNAtools
            [0]
                1) Find DNA sequence in the original sequence.
                2) Compare whether the first or second amino acid sequence is: greater, lesser, equal, greater than or equal, lesser than or equal.
                3) Know the percentage that two DNA sequences have in common.
                4) Generate permutations from the DNA sequence.
                5) Generate a being from bits.

        """)

        option2 = input(str("Please select a option: "))

        if option2 == "1":
        #· Find DNA sequence in the original sequence.
            str1 = input("Input the DNA chain (A, C, G, T): ")
            characters_str1 = str(len(str1))
            print("This sequence have " + characters_str1 + " characters")

            str2 = input("Insert the sequence you want to find: ")
            characters_str2 = str(len(str2))
            print("This sequence have " + characters_str2 + " characters")

            str1 = str1.replace(str2,f'\033[96m{str2}\033[0m')
            print(str1)

        elif option2 == "2":
        #· Compare whether the first or second amino acid sequence is: greater, lesser, equal, greater than or equal, lesser than or equal.
            DNA1 = input("Input the first DNA chain (A, C, G, T): ")
            characters_str1 = str(len(DNA1))
            print("This sequence have " + characters_str1 + " characters")
            DNA2 = input("Input the second DNA chain (A, C, G, T): ")
            characters_str2 = str(len(DNA2))
            print("This sequence have " + characters_str2 + " characters")

            print('The first DNA sequence is greater than the second DNA sequence:', DNA1.upper() > DNA2.upper())
            print('The first DNA sequence is lesser than the second DNA sequence:', DNA1.upper() < DNA2.upper())
            print('The first DNA sequence is equal to the second DNA sequence:', DNA1.upper() == DNA2.upper())
            print('The first DNA sequence is greater than or equal to the second DNA sequence:', DNA1.upper() >= DNA2.upper())
            print('The first DNA sequence is lesser than or equal to the second DNA sequence:', DNA1.upper() <= DNA2.upper())
        
        elif option2 == "3":
        #· This function returns the match percentage of two DNA sequences.

            first = str(input("Input the dna sequence (A, G, C, T): "))
            characters_str1 = str(len(first))
            print("This sequence have " + characters_str1 + " characters")
            second = str(input("Input other dna sequence to campare with another sequence: "))
            characters_str2 = str(len(second))
            print("This sequence have " + characters_str2 + " characters")

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
            DNA = input("Input the DNA sequence (A, G, C, T): ")
            characters_str = str(len(DNA))
            print("This sequence have " + characters_str + " characters")
            permutation = int(input("By which number do you want to permutation: "))

            letras = [letra for letra in DNA if letra in 'ACTG' or 'actg']
            print(letras)
            print()
            permutaciones = [p for p in permutations(letras, permutation)]
            contador = 0

            for p in permutaciones:
                print(p)
                contador += 1
            print()
            print(DNA)
            print(f'Number of permutations generated: {contador}.')
        elif option2 == "5":
            
            index_nitrogenous_bases = {"Adenina": 0, "Citosina": 1, "Guanina": 2, "Timina": 3}
                
            nitrogenous_bases = [
                #Quarks
                {'Name': 'Adenina', 'Initial letter' : 'A', 'Binary': "00000000", 'Color_name' : '\033[33mAdenina\033[0m', 'Color_binary' : '\033[33m00000000\033[0m'},
                {'Name': 'Citosina', 'Initial letter' : 'C', 'Binary': "00000001", 'Color_name' : '\033[34mCitosina\033[0m', 'Color_binary' : '\033[34m00000001\033[0m'},
                {'Name': 'Guanina', 'Initial letter' : 'G', 'Binary': "00000010", 'Color_name' : '\033[32mGuanina\033[0m', 'Color_binary' : '\033[32m00000010\033[0m'},
                {'Name': 'Timina', 'Initial letter' : 'T', 'Binary': "00000011", 'Color_name' : '\033[31mTimina\033[0m', 'Color_binary' : '\033[31m00000011\033[0m'},]


            bits = ['0','1']

            n = int(input("Enter the number of bits you want to generate: "))
            total_bytes = n/8

            # We create a string that returns zeros and ones randomly.
            def random_binary(n):
                return ''.join(random.choice(bits) for _ in range(n))

            # Every eight bits is a byte and we separate it with a colon (:).
            def binary_to_bytes(binary):
                return ':'.join(binary[i:i+8] for i in range(0, len(binary), 8))

            # We store the bits in a string
            bytes = binary_to_bytes(random_binary(n))

            print("\nI show you the bytes that have coincided with some nitrogenous bases:\n")

            for byte in bytes.split(':'):
                for nitrogen_base in nitrogenous_bases:
                    if nitrogen_base['Binary'] == byte:
                        # We change the byte for the colored byte.
                        bytes = bytes.replace(byte, nitrogen_base['Color_binary'])
                        break
            print(bytes)

            print("\nI show you the name of the nitrogenous bases that has matched:\n")
            for byte_colored in range(len(bytes.split(':'))):
                for nitrogen_base in nitrogenous_bases:
                    if nitrogen_base['Color_binary'] in bytes.split(':')[byte_colored]:
                        # We change the colored byte for the colored name.
                        bytes = bytes.replace(bytes.split(':')[byte_colored], nitrogen_base['Color_name'])
                        break
            print(bytes)
            print("\nI have created", total_bytes, "bytes.")
        elif option2 == "0":
            print("Until next time.")
            exit()
        else:
            print("Invalid option.")
            print("Please, try again.")
    elif option == "0":
        print("Until next time.")
        exit()
    else:
        print("Invalid option.")
        print("Please, try again.")
