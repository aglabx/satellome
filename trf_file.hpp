#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct TR {
    std::string tr_id;
    std::string tr_array;
    std::vector<size_t> embeding;
    std::vector<double> norm_embeding;
};


size_t read_trs_file_to_vector(const std::string &trs_file,
                            std::vector<TR> &trs_vector) {
    std::ifstream trs_stream(trs_file);

    if (!trs_stream.is_open()) {
        std::cerr << "Can't open file: " << trs_file << std::endl;
        exit(1);
    }


    std::string str;
    while(std::getline(trs_stream, str)) {
        
        std::stringstream buffer(str);
        std::string temp;
        std::vector<std::string> values;

        while(std::getline(buffer, temp, '\t') ) {
            values.push_back(temp);
        }

        if (values.size() < 15) {
            continue;
        }

        TR tr;
        tr.tr_id = values[1];
        tr.tr_array = values[14];
        trs_vector.push_back(tr);
    }

    return trs_vector.size();
}