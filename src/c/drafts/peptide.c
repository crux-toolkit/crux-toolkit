typedef struct Peptide
{
    char PrefixAmino; // The base BEFORE the peptide starts.  (Useful for checking trypsin fragmentation)
    char SuffixAmino; // The base AFTER the peptide starts.
    int AminoIndex[MAX_PT_MODS]; 
    MassDelta* ModType[MAX_PT_MODS];
    struct PeptideMatch* First;
} Peptide;

