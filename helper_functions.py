import numpy as np

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    str1 = reference_sequence
    str2 = distant_sequence
    if len(str1) != len(str2):
        raise ValueError("Input strings must have the same length")

    # Create filtered strings without positions with indels ("-")
    filtered_str1 = ''.join(c1 if c2 != '-' else '-' for c1, c2 in zip(str1, str2))
    filtered_str1 = filtered_str1.replace("-","")
    filtered_str2 = ''.join(c2 if c1 != '-' else '-' for c1, c2 in zip(str1, str2))
    filtered_str2 = filtered_str2.replace("-","")

    #print(filtered_str1)
    #print(filtered_str2)

    # Calculate Hamming distance on the filtered strings
    distance = sum(c1 != c2 for c1, c2 in zip(filtered_str1, filtered_str2))
    #print(distance)
    
    p = distance / len(filtered_str1)
    #print(p)
    corrected_distance = -3/4 * np.log(1 - (4/3) * p)

    return corrected_distance*len(filtered_str1)


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    referene_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()