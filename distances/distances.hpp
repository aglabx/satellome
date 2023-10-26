#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

struct Distance {
    size_t tr_idA;
    size_t tr_idB;
    double distance;
};

double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0;
    for (size_t i = 0; i < v1.size(); i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

double norm(const std::vector<double>& v) {
    double sum = 0;
    for (size_t i = 0; i < v.size(); i++) {
        sum += v[i] * v[i];
    }
    return std::sqrt(sum);
}

double cosine_distance(const std::vector<double>& v1, const std::vector<double>& v2) {
    double dotProd = dot_product(v1, v2);
    double norm1 = norm(v1);
    double norm2 = norm(v2);
    double cosineSim = dotProd / (norm1 * norm2);
    return 1 - cosineSim;
}