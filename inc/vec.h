
template <typename T>
void operator*=(std::vector<T> &p, const double v) {
    for (std::size_t i = 0; i < p.size(); ++i) {
        p[i] *= v;
    }
}

template<typename T>
inline std::vector<T> operator*(const T v, const std::vector<T>& q) {
    // Create a new vector to store the result
    std::vector<T> result(q.size());

    // Add the corresponding elements of p and q and store the result in the result vector
    for (std::size_t i = 0; i < q.size(); ++i) {
        result[i] = v * q[i];
    }

    return result;
}

template<typename T>
inline std::vector<T> operator*(const std::vector<T>& q, const T v) {
    return v*q;
}

template<typename T>
inline std::vector<T> operator+(const std::vector<T>& p, const std::vector<T>& q) {
    // Make sure both vectors have the same size
    if (p.size() != q.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }

    // Create a new vector to store the result
    std::vector<T> result(p.size());

    // Add the corresponding elements of p and q and store the result in the result vector
    for (std::size_t i = 0; i < p.size(); ++i) {
        result[i] = p[i] + q[i];
    }

    return result;
}

template<typename T>
inline std::vector<T> operator-(const std::vector<T>& p, const std::vector<T>& q) {
    // Make sure both vectors have the same size
    if (p.size() != q.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }

    // Create a new vector to store the result
    std::vector<T> result(p.size());

    // Add the corresponding elements of p and q and store the result in the result vector
    for (std::size_t i = 0; i < p.size(); ++i) {
        result[i] = p[i] - q[i];
    }

    return result;
}