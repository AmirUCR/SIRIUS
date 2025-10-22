#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "logging.hpp"
#include "tables.hpp"

#include "csv.hpp"

namespace sirius
{
    SIRIUSTables::SIRIUSTables()
    : skip_aa{}
    , invariant_codon_table{}
    , reduced_codon_lengths_table{}
    , reduced_codon_table{}
    , translate_codon_table{}
    {
        skip_aa = {};
        
        invariant_codon_table = {
            {'A', "___"}, // {"GCT", "GCC", "GCA", "GCG"}
            {'C', "___"}, // {"TGT", "TGC"}
            {'D', "___"}, // {"GAT", "GAC"}
            {'E', "___"}, // {"GAA", "GAG"}
            {'F', "___"}, // {"TTT", "TTC"}
            {'G', "___"}, // {"GGT", "GGC", "GGA", "GGG"}
            {'H', "___"}, // {"CAT", "CAC"}
            {'I', "___"}, // {"ATT", "ATC", "ATA"}
            {'K', "___"}, // {"AAA", "AAG"}
            {'L', "___"}, // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
            {'M', "___"},
            {'N', "___"}, // {"AAT", "AAC"}
            {'P', "___"}, // {"CCT", "CCC", "CCA", "CCG"}
            {'Q', "___"}, // {"CAA", "CAG"}
            {'R', "___"}, // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
            {'S', "___"}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
            {'T', "___"}, // {"ACT", "ACC", "ACA", "ACG"}
            {'V', "___"}, // {"GTT", "GTC", "GTA", "GTG"}
            {'W', "___"}, // TGG
            {'Y', "___"}, // {"TAT", "TAC"}
            {'*', "___"}  // {"TAA", "TAG", "TGA"}
        };

        reduced_codon_lengths_table = {
            {'A', 3},
            {'C', 3},
            {'D', 3},
            {'E', 3},
            {'F', 3},
            {'G', 3},
            {'H', 3},
            {'I', 3},
            {'K', 3},
            {'L', 3},
            {'M', 3},
            {'N', 3},
            {'P', 3},
            {'Q', 3},
            {'R', 3},
            {'S', 3},
            {'T', 3},
            {'V', 3},
            {'W', 3},
            {'Y', 3},
            {'*', 3}
        };

        reduced_codon_table = {
            {'A', {"GCA", "GCC", "GCG", "GCT"}},
            {'C', {"TGC", "TGT"}},
            {'D', {"GAC", "GAT"}},
            {'E', {"GAA", "GAG"}},
            {'F', {"TTC", "TTT"}},
            {'G', {"GGA", "GGC", "GGG", "GGT"}},
            {'H', {"CAC", "CAT"}},
            {'I', {"ATA", "ATC", "ATT"}},
            {'K', {"AAA", "AAG"}},
            {'L', {"CTA", "CTC", "CTG", "CTT", "TTA", "TTG"}},
            {'M', {"ATG"}},
            {'N', {"AAC", "AAT"}},
            {'P', {"CCA", "CCC", "CCG", "CCT"}},
            {'Q', {"CAA", "CAG"}},
            {'R', {"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"}},
            {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}},
            {'T', {"ACA", "ACC", "ACG", "ACT"}},
            {'V', {"GTA", "GTC", "GTG", "GTT"}},
            {'W', {"TGG"}},
            {'Y', {"TAC", "TAT"}},
            {'*', {"TAA", "TAG", "TGA"}}
        };

        translate_codon_table = {
            // Phenylalanine
            {"TTT", 'F'},
            {"TTC", 'F'},
            // Leucine
            {"TTA", 'L'},
            {"TTG", 'L'},
            {"CTT", 'L'},
            {"CTC", 'L'},
            {"CTA", 'L'},
            {"CTG", 'L'},
            // Isoleucine
            {"ATT", 'I'},
            {"ATC", 'I'},
            {"ATA", 'I'},
            // Methionine (Start)
            {"ATG", 'M'},
            // Valine
            {"GTT", 'V'},
            {"GTC", 'V'},
            {"GTA", 'V'},
            {"GTG", 'V'},
            // Serine
            {"TCT", 'S'},
            {"TCC", 'S'},
            {"TCA", 'S'},
            {"TCG", 'S'},
            {"AGT", 'S'},
            {"AGC", 'S'},
            // Proline
            {"CCT", 'P'},
            {"CCC", 'P'},
            {"CCA", 'P'},
            {"CCG", 'P'},
            // Threonine
            {"ACT", 'T'},
            {"ACC", 'T'},
            {"ACA", 'T'},
            {"ACG", 'T'},
            // Alanine
            {"GCT", 'A'},
            {"GCC", 'A'},
            {"GCA", 'A'},
            {"GCG", 'A'},
            // Tyrosine
            {"TAT", 'Y'},
            {"TAC", 'Y'},
            // Histidine
            {"CAT", 'H'},
            {"CAC", 'H'},
            // Glutamine
            {"CAA", 'Q'},
            {"CAG", 'Q'},
            // Asparagine
            {"AAT", 'N'},
            {"AAC", 'N'},
            // Lysine
            {"AAA", 'K'},
            {"AAG", 'K'},
            // Aspartic Acid
            {"GAT", 'D'},
            {"GAC", 'D'},
            // Glutamic Acid
            {"GAA", 'E'},
            {"GAG", 'E'},
            // Cysteine
            {"TGT", 'C'},
            {"TGC", 'C'},
            // Tryptophan
            {"TGG", 'W'},
            // Arginine
            {"CGT", 'R'},
            {"CGC", 'R'},
            {"CGA", 'R'},
            {"CGG", 'R'},
            {"AGA", 'R'},
            {"AGG", 'R'},
            // Glycine
            {"GGT", 'G'},
            {"GGC", 'G'},
            {"GGA", 'G'},
            {"GGG", 'G'},
            // Stop codons
            {"TAA", '*'},
            {"TAG", '*'},
            {"TGA", '*'}
        };
    };

    void SIRIUSTables::build_rscu_map_from_csv(const std::string &filename)
    {
        csv::CSVReader reader(filename);

        double rscu;
        std::string codon;
        std::string aa_str;

        // The 'reader' yields 'csv::CSVRow' objects.
        for (csv::CSVRow& row : reader)
        {
            try 
            {
                // Access columns by name and retrieve the value, converting to the required type.
                aa_str = row["AmOneLet"].get<std::string>();
                codon = row["Codon"].get<std::string>();
                rscu = row["RSCU"].get<double>();
            } 
            catch (const std::exception& e) 
            {
                // Catch errors during type conversion or if a required column is missing.
                logg.throw_formatted_error("Error: Malformed data in " + filename + ". Detail: " + e.what());
            }

            if (aa_str.empty() || codon.size() != 3)
            {
                logg.throw_formatted_error("Error: Malformed CSV line in " + filename);
            }

            char aa = aa_str[0];
            auto it = this->invariant_codon_table.find(aa);
            if (it == this->invariant_codon_table.end())
            {
                logg.throw_formatted_error("Error: Unknown amino acid '" + std::string(1, aa) + "'.");
            }

            const std::string &pattern = it->second;
            std::string varying_part;
            for (size_t i = 0; i < pattern.size(); ++i)
            {
                if (pattern[i] == '_')
                {
                    varying_part += codon[i];
                }
            }

            if (!varying_part.empty())
            {
                this->rscu_map[aa][varying_part] = rscu;
            }
        }
    }
}
