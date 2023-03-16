// g++ -std=c++2a -pthread -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -O3 -rdynamic json.hpp bpe.cpp -o bpe.exe 

#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>
#include <numeric>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include "json.hpp"


using json = nlohmann::json;

std::unordered_map<std::string, int> alphabet = {
    {"A", 1},
    {"C", 2},
    {"G", 3},
    {"T", 4},
    {"N", 5},
    {"~", 0}
};

// Custom hash function for std::tuple
template<typename T>
struct tuple_hash {
    template<typename Tuple, size_t N>
    struct hasher {
        static void hash(Tuple const& t, std::size_t& seed) {
            hasher<Tuple, N - 1>::hash(t, seed);
            seed ^= std::hash<typename std::tuple_element<N - 1, Tuple>::type>()(std::get<N - 1>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    };

    template<typename Tuple>
    struct hasher<Tuple, 1> {
        static void hash(Tuple const& t, std::size_t& seed) {
            seed ^= std::hash<typename std::tuple_element<0, Tuple>::type>()(std::get<0>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    };

    template<typename... Args>
    std::size_t operator()(std::tuple<Args...> const& t) const {
        std::size_t seed = 0;
        hasher<std::tuple<Args...>, sizeof...(Args)>::hash(t, seed);
        return seed;
    }
};

typedef std::tuple<size_t, size_t> kmer;

void get_sequences(std::string& file_name, std::vector<std::string>& seqs) {
    std::ifstream fh(file_name);
    if (fh.is_open()) {
        std::string line;
        while (std::getline(fh, line)) {
            std::vector<std::string> d;
            size_t pos = 0;
            while ((pos = line.find('\t')) != std::string::npos) {
                d.push_back(line.substr(0, pos));
                line.erase(0, pos + 1);
            }
            seqs.push_back(d[14]);
        }
        fh.close();
    }
}

std::string get_dataset(std::vector<std::string>& seqs) {
    std::stringstream ss;
    bool first = true;
    for (const auto& s : seqs) {
        if (!first) {
            ss << "~";
        } else {
            first = false;
        }
        ss << s;
    }
    return ss.str();
}

std::vector<size_t> convert_to_vector(std::string& dataset) {
    std::vector<size_t> seq;
    seq.reserve(dataset.size()); // Reserve space
    for (const auto& x : dataset) {
        seq.push_back(alphabet[std::string(1, x)]);
    }
    return seq;
}

// std::string tokens_to_letters(const std::unordered_map<std::string, int>& alphabet, const std::unordered_map<size_t, kmer>& tokens, const std::vector<size_t>& seq, bool recursed = false) {
//     std::unordered_map<int, std::string> reverse_alphabet;
//     for (const auto& entry : alphabet) {
//         reverse_alphabet[entry.second] = entry.first;
//     }

//     std::stringstream ss;
//     for (const auto& token : seq) {
//         if (tokens.find(token) != tokens.end()) {
//             // This is a merged token (kmer)
//             kmer kmer_ = tokens.at(token);
//             if (!recursed) {
//                 // Recursively handle merged tokens
//                 ss << tokens_to_letters(alphabet, tokens, kmer_, true);
//             } else {
//                 for (const auto& element : kmer_) {
//                     ss << reverse_alphabet[element];
//                 }
//             }
//         } else {
//             // This is a single letter from the alphabet
//             ss << reverse_alphabet[token];
//         }
//     }

//     return ss.str();
// }

// void output_json(const std::unordered_map<std::string, int>& alphabet, const std::unordered_map<size_t, kmer>& tokens, const std::vector<size_t>& seq, const std::string& output_file) {
//     json output;

//     output["version"] = "1.0";
//     output["truncation"] = nullptr;
//     output["padding"] = nullptr;
//     // Add other fields, such as "added_tokens", "normalizer", "pre_tokenizer", "post_processor", "decoder" and "model.type"

//     // Populate the "vocab" object
//     json vocab;
//     for (const auto& entry : alphabet) {
//         vocab[entry.first] = entry.second;
//     }

//     // Populate the "merges" array
//     json merges = json::array();
//     for (const auto& entry : tokens) {
//         std::stringstream ss;
//         for (size_t i = 0; i < entry.second.size(); i++) {
//             ss << alphabet.at(entry.second[i]);
//             if (i < entry.second.size() - 1) {
//                 ss << " ";
//             }
//         }
//         merges.push_back(ss.str());
//     }

//     output["model"]["vocab"] = vocab;
//     output["model"]["merges"] = merges;

//     // Write the output to a file
//     std::ofstream out_file(output_file);
//     if (out_file.is_open()) {
//         out_file << output.dump(2);  // Use '2' for pretty-printing with a 2-space indent
//         out_file.close();
//     }
// }


int main() {
    std::string file_name = "/mnt/data/podgornaya/rana_temporaria/users/akomissarov/trf/GCF_905171775.1_aRanTem1.1_genomic.1kb.trf";
    std::string output_file = "/mnt/data/podgornaya/rana_temporaria/users/akomissarov/trf/GCF_905171775.1_aRanTem1.1_genomic.1kb.bpe.json";
    std::vector<std::string> seqs;
    std::cout << "read file" << std::endl;
    get_sequences(file_name, seqs);
    std::cout << "get dataset" << std::endl;
    std::string dataset = get_dataset(seqs);
    
    std::vector<size_t> seq = convert_to_vector(dataset);

    std::vector<kmer> merged;
    int k = 2;
    int L = alphabet.size() + 1;
    std::unordered_map<size_t, kmer> tokens;
    std::unordered_map<kmer, size_t, tuple_hash<size_t>> rev_tokens;

    std::vector<size_t> new_seq;
    std::vector<bool> to_replace(seq.size(), false);
    while (true) {
        std::cout << "Tokens " << L << " count reps ";

        std::vector<std::vector<size_t>> c(L, std::vector<size_t>(L, 0));
        for (size_t i = 0; i < seq.size() - k + 1; i++) {
            if (seq[i] && seq[i+1]) {
                c[seq[i]][seq[i+1]]++;
            }
        }

        std::cout << "find max ";

        size_t max_count = 0;
        kmer rep;
        for (size_t i = 1; i < c.size(); ++i) {
            for (size_t j = 1; j < c[i].size(); ++j) {
                if (c[i][j] > max_count) {
                    max_count = c[i][j];
                    rep = std::make_tuple(i, j);
                }
            }
        }
        size_t tf = max_count;
        if (tf == 1) {
            break;
        }
        // if (L > 100) {
        //     break;
        // }

        std::cout << std::get<0>(rep) << " " << std::get<1>(rep) << " " << tf << " : ";
        merged.push_back(rep);
        tokens[L] = rep;
        rev_tokens[rep] = L;
        
        std::fill(to_replace.begin(), to_replace.end(), false);
        for (size_t i = 0; i < seq.size() - k + 1; i++) {
            if (seq[i] && seq[i+1]) {
                kmer kmer_ = std::make_tuple(seq[i], seq[i+1]);
                if (kmer_ == rep) {
                    to_replace[i] = true;
                }
            }
        }
        
        std::cout << "replace: ";
        new_seq.clear();
        new_seq.reserve(seq.size());
        for (size_t i = 0; i < seq.size();) {
            if (to_replace[i]) {
                new_seq.push_back(L);
                i += 2;
            } else {
                new_seq.push_back(seq[i]);
                i += 1;
            }
        }
        L += 1;
        std::cout << seq.size() << " -> " << new_seq.size();
        std::cout << " new seq copy ";
        seq = std::vector<size_t>(new_seq.begin(), new_seq.end());
        std::cout << "done" << std::endl;
    }

    std::ofstream out_file(output_file);
    if (out_file.is_open()) {
        out_file << "alphaber" << std::endl;
        for (const auto& element : alphabet) {
            out_file << element.first << " " << element.second << std::endl;
        }
        out_file << "merged" << std::endl;
        for (const auto& element : merged) {
            out_file << std::get<0>(element) << " " << std::get<1>(element) << std::endl;
        }
        out_file << "tokens" << std::endl;
        for (const auto& element : tokens) {
            out_file << element.first << " " << std::get<0>(element.second) << " " << std::get<1>(element.second) << std::endl;
        }
        out_file << "sequence" << std::endl;
        for (const auto& element : seq) {
            out_file << element << " ";
        }
        out_file << std::endl;
        out_file.close();
    }


    return 0;
}
